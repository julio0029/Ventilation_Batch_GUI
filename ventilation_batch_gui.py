"""Local batch ventilation analysis for shark and ray videos.

Select a locally synced Dropbox folder. The application scans all video files,
estimates periodic motion in a normalised ROI, and writes one CSV row per file.
"""

from __future__ import annotations

import csv
import queue
import threading
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Callable

import cv2
import numpy as np
import tkinter as tk
from tkinter import filedialog, messagebox, ttk

try:
    import neurokit2 as nk
except ImportError:  # Message below explains how to install the required dependency.
    nk = None


VIDEO_EXTENSIONS = {".mp4", ".mov", ".avi", ".mkv", ".m4v", ".wmv"}
DEFAULT_ROI = (0.50, 1.00, 0.15, 0.40)  # y start, y end, x start, x end
SPECIES = ("epaulette_shark", "shovelnose_ray")


@dataclass
class Settings:
    rois: dict[str, tuple[float, float, float, float]]
    min_bpm: float
    max_bpm: float
    prominence: float
    diagnostics: bool


def safe_metadata(relative_video: Path, root_name: str) -> dict[str, str]:
    """Infer experimental labels from any folder order using case insensitive keywords."""
    tokens = [part.lower() for part in relative_video.parts[:-1]]
    text = " ".join(tokens + [relative_video.stem.lower()])

    def label(mapping: dict[str, str]) -> str:
        for key, value in mapping.items():
            if key in text:
                return value
        return "unclassified"

    if any(term in text for term in ("epaulette", "hemiscyllium", "ocellatum")):
        species = "epaulette_shark"
        species_group = "shark"
    elif any(term in text for term in ("shovelnose", "aptychotrema", "rostrata")):
        species = "shovelnose_ray"
        species_group = "ray"
    elif "shark" in text:
        species = "unclassified_shark"
        species_group = "shark"
    elif "ray" in text:
        species = "unclassified_ray"
        species_group = "ray"
    else:
        species = "unclassified"
        species_group = "unclassified"

    parent = relative_video.parent
    return {
        "root_folder": root_name,
        "top_level_folder": relative_video.parts[0] if relative_video.parts else "",
        "relative_folder": str(parent) if str(parent) != "." else "",
        "immediate_folder": parent.name if str(parent) != "." else "",
        "file_name": relative_video.name,
        "species": species,
        "species_group": species_group,
        "condition": label({"hypoxia": "hypoxia", "normoxia": "normoxia", "normaxia": "normoxia"}),
        "treatment": label({"pbs": "PBS", "poly": "Poly"}),
    }


def extract_motion(video_path: Path, metadata: dict[str, str], settings: Settings) -> tuple[np.ndarray, np.ndarray, float, int, int]:
    cap = cv2.VideoCapture(str(video_path))
    if not cap.isOpened():
        raise RuntimeError("OpenCV could not open this video")
    fps = float(cap.get(cv2.CAP_PROP_FPS))
    width = int(cap.get(cv2.CAP_PROP_FRAME_WIDTH))
    height = int(cap.get(cv2.CAP_PROP_FRAME_HEIGHT))
    frames = int(cap.get(cv2.CAP_PROP_FRAME_COUNT))
    if fps <= 0 or width <= 0 or height <= 0 or frames < 3:
        cap.release()
        raise RuntimeError("Video metadata is incomplete or the video is too short")

    species = metadata["species"]
    if species not in settings.rois:
        cap.release()
        raise RuntimeError("Species is not recognised. Name folders or files epaulette_shark or shovelnose_ray before analysis")
    y0f, y1f, x0f, x1f = settings.rois[species]
    y0, y1 = int(y0f * height), int(y1f * height)
    x0, x1 = int(x0f * width), int(x1f * width)
    if y1 <= y0 or x1 <= x0:
        cap.release()
        raise RuntimeError("ROI has zero area")

    previous = None
    motion, timestamps = [], []
    index = 0
    while True:
        ok, frame = cap.read()
        if not ok:
            break
        gray = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)
        roi = gray[y0:y1, x0:x1].astype(np.float32)
        if previous is None:
            motion.append(0.0)
        else:
            motion.append(float(np.mean(np.abs(roi - previous))))
        timestamps.append(index / fps)
        previous = roi
        index += 1
    cap.release()
    return np.asarray(timestamps), np.asarray(motion), fps, width, height


def moving_average(values: np.ndarray, window: int) -> np.ndarray:
    """Centered moving average using NumPy only."""
    if window <= 1:
        return values.copy()
    pad_left = window // 2
    pad_right = window - 1 - pad_left
    padded = np.pad(values, (pad_left, pad_right), mode="edge")
    return np.convolve(padded, np.ones(window) / window, mode="valid")


def linear_detrend(values: np.ndarray) -> np.ndarray:
    positions = np.arange(len(values), dtype=float)
    slope, intercept = np.polyfit(positions, values, 1)
    return values - (slope * positions + intercept)


def neurokit_peaks(values: np.ndarray, distance: int, threshold: float) -> np.ndarray:
    """Use NeuroKit2 to identify candidate contractions, then enforce cycle spacing."""
    if nk is None:
        raise RuntimeError("NeuroKit2 is required. Run python -m pip install -r requirements.txt in a clean Python environment")
    candidates = np.asarray(
        nk.signal_findpeaks(values, relative_height_min=threshold, relative_mean=True)["Peaks"], dtype=int
    )
    if len(candidates) == 0:
        return np.asarray([], dtype=int)
    selected: list[int] = []
    for index in sorted(candidates, key=lambda item: values[item], reverse=True):
        if all(abs(index - existing) >= distance for existing in selected):
            selected.append(index)
    return np.asarray(sorted(selected), dtype=int)


def estimate_rate(timestamps: np.ndarray, motion: np.ndarray, fps: float, settings: Settings) -> dict[str, float | int | str]:
    if len(motion) < max(30, int(fps * 5)):
        raise RuntimeError("Video is too short for a stable ventilation estimate")
    smooth_window = max(2, int(round(fps * 0.07)))
    smoothed = moving_average(motion, smooth_window)
    detrended = linear_detrend(smoothed)
    standard_deviation = float(np.std(detrended))
    if standard_deviation == 0:
        raise RuntimeError("ROI motion signal is constant")
    scaled = detrended / standard_deviation
    window = np.hanning(len(scaled))
    freqs = np.fft.rfftfreq(len(scaled), d=1 / fps)
    psd = np.abs(np.fft.rfft(scaled * window)) ** 2
    band = (freqs >= settings.min_bpm / 60) & (freqs <= settings.max_bpm / 60)
    if not np.any(band):
        raise RuntimeError("Video sampling rate cannot resolve the requested ventilation range")
    band_freqs, band_psd = freqs[band], psd[band]
    dominant_hz = float(band_freqs[np.argmax(band_psd)])
    fft_bpm = dominant_hz * 60
    distance = max(1, int(0.70 * fps / dominant_hz)) if dominant_hz > 0 else max(1, int(fps))
    peaks = neurokit_peaks(scaled, distance=distance, threshold=settings.prominence)
    if len(peaks) >= 2:
        intervals = np.diff(timestamps[peaks])
        mean_interval = float(np.mean(intervals))
        peak_bpm = 60 / mean_interval if mean_interval > 0 else float("nan")
        interval_cv = float(np.std(intervals) / mean_interval) if mean_interval > 0 else float("nan")
        final_bpm = float(np.mean([fft_bpm, peak_bpm]))
        method = "FFT and NeuroKit2 peak consensus"
    else:
        peak_bpm = float("nan")
        interval_cv = float("nan")
        final_bpm = fft_bpm
        method = "FFT only; insufficient NeuroKit2 peaks"
    return {
        "ventilation_bpm": final_bpm,
        "fft_bpm": fft_bpm,
        "peak_bpm": peak_bpm,
        "peak_count": int(len(peaks)),
        "inter_peak_cv": interval_cv,
        "method": method,
        "signal_sd": standard_deviation,
        "motion": motion,
        "smoothed": smoothed,
        "scaled": scaled,
        "peaks": peaks,
    }


def save_diagnostic(video_path: Path, output_path: Path, timestamps: np.ndarray, result: dict, metadata: dict) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure, axes = plt.subplots(2, 1, figsize=(12, 7), sharex=True)
    axes[0].plot(timestamps, result["motion"], color="#b0bec5", linewidth=0.6, label="Raw motion")
    axes[0].plot(timestamps, result["smoothed"], color="#1565c0", linewidth=1, label="Smoothed motion")
    axes[0].set_ylabel("Mean frame difference")
    axes[0].legend()
    axes[1].plot(timestamps, result["scaled"], color="#37474f", linewidth=0.7, label="Detrended motion")
    peaks = result["peaks"]
    if len(peaks):
        axes[1].plot(timestamps[peaks], result["scaled"][peaks], "v", color="#d32f2f", markersize=5, label="Detected breaths")
    axes[1].set_ylabel("Normalised signal")
    axes[1].set_xlabel("Time, seconds")
    axes[1].legend()
    figure.suptitle(f"{metadata['file_name']} | {result['ventilation_bpm']:.1f} breaths per minute")
    figure.tight_layout()
    figure.savefig(output_path, dpi=150)
    plt.close(figure)


def analyse_video(video_path: Path, root: Path, settings: Settings, diagnostics_dir: Path | None) -> dict[str, object]:
    relative = video_path.relative_to(root)
    metadata = safe_metadata(relative, root.name)
    try:
        timestamps, motion, fps, width, height = extract_motion(video_path, metadata, settings)
        result = estimate_rate(timestamps, motion, fps, settings)
        if diagnostics_dir is not None:
            diagnostics_dir.mkdir(parents=True, exist_ok=True)
            diagnostic_name = "__".join(relative.with_suffix("").parts) + ".png"
            save_diagnostic(video_path, diagnostics_dir / diagnostic_name, timestamps, result, metadata)
        return {
            **metadata,
            "video_path": str(video_path),
            "fps": round(fps, 4),
            "frame_width": width,
            "frame_height": height,
            "duration_seconds": round(float(timestamps[-1]), 3),
            "ventilation_bpm": round(float(result["ventilation_bpm"]), 3),
            "fft_bpm": round(float(result["fft_bpm"]), 3),
            "peak_bpm": "" if np.isnan(result["peak_bpm"]) else round(float(result["peak_bpm"]), 3),
            "peak_count": result["peak_count"],
            "inter_peak_cv": "" if np.isnan(result["inter_peak_cv"]) else round(float(result["inter_peak_cv"]), 4),
            "peak_algorithm": "NeuroKit2 signal_findpeaks",
            "method": result["method"],
            "signal_sd": round(float(result["signal_sd"]), 5),
            "status": "ok",
            "error": "",
        }
    except Exception as exc:
        return {**metadata, "video_path": str(video_path), "status": "failed", "error": str(exc)}


class VentilationGUI(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Shark and Ray Ventilation Batch Analysis")
        self.geometry("850x680")
        self.queue: queue.Queue = queue.Queue()
        self.root_var = tk.StringVar()
        self.output_var = tk.StringVar()
        self.roi_vars = {
            species: [tk.StringVar(value=f"{x:.2f}") for x in DEFAULT_ROI]
            for species in SPECIES
        }
        self.min_bpm_var = tk.StringVar(value="20")
        self.max_bpm_var = tk.StringVar(value="180")
        self.prominence_var = tk.StringVar(value="0.0")
        self.diagnostics_var = tk.BooleanVar(value=True)
        self._build()
        self.after(150, self._drain_queue)

    def _build(self):
        frame = ttk.Frame(self, padding=14)
        frame.pack(fill="both", expand=True)
        ttk.Label(frame, text="Shark and Ray Ventilation Batch Analysis", font=("Arial", 16, "bold")).grid(row=0, column=0, columnspan=4, sticky="w")
        ttk.Label(frame, text="Select the locally synced Dropbox folder that contains your video subfolders.").grid(row=1, column=0, columnspan=4, sticky="w", pady=(4, 12))
        self._path_row(frame, 2, "Video root folder", self.root_var, self.choose_root)
        self._path_row(frame, 3, "Results CSV", self.output_var, self.choose_output, save=True)
        ttk.Separator(frame).grid(row=4, column=0, columnspan=4, sticky="ew", pady=12)
        ttk.Label(frame, text="Species specific ROI fractions: y start, y end, x start, x end").grid(row=5, column=0, columnspan=5, sticky="w")
        for row, species in enumerate(SPECIES, 6):
            ttk.Label(frame, text=species.replace("_", " ").title()).grid(row=row, column=0, sticky="w")
            for index, variable in enumerate(self.roi_vars[species]):
                ttk.Entry(frame, textvariable=variable, width=8).grid(row=row, column=index + 1, sticky="w", padx=3)
            ttk.Button(frame, text="Choose ROI", command=lambda value=species: self.choose_roi(value)).grid(row=row, column=5, sticky="w", padx=(6, 0))
        ttk.Label(frame, text="Select a separate ROI for each species. Use comparable camera framing within each species.").grid(row=8, column=0, columnspan=6, sticky="w", pady=(6, 0))
        ttk.Label(frame, text="Ventilation range, breaths per minute").grid(row=9, column=0, sticky="w", pady=(10, 0))
        ttk.Entry(frame, textvariable=self.min_bpm_var, width=8).grid(row=9, column=1, sticky="w", padx=3, pady=(10, 0))
        ttk.Entry(frame, textvariable=self.max_bpm_var, width=8).grid(row=9, column=2, sticky="w", padx=3, pady=(10, 0))
        ttk.Label(frame, text="NeuroKit2 relative peak height threshold").grid(row=10, column=0, columnspan=2, sticky="w", pady=(6, 0))
        ttk.Entry(frame, textvariable=self.prominence_var, width=8).grid(row=10, column=2, sticky="w", padx=3, pady=(6, 0))
        ttk.Checkbutton(frame, text="Save one diagnostic plot per video", variable=self.diagnostics_var).grid(row=11, column=0, columnspan=3, sticky="w", pady=8)
        ttk.Button(frame, text="Analyse all videos", command=self.start).grid(row=12, column=0, sticky="w", pady=8)
        self.progress = ttk.Progressbar(frame, mode="determinate")
        self.progress.grid(row=12, column=1, columnspan=5, sticky="ew", padx=8)
        ttk.Label(frame, text="Metadata require epaulette_shark or shovelnose_ray, hypoxia or normoxia, and PBS or Poly in folders or filenames.").grid(row=13, column=0, columnspan=6, sticky="w", pady=(6, 2))
        self.log = tk.Text(frame, height=21, wrap="word")
        self.log.grid(row=14, column=0, columnspan=6, sticky="nsew", pady=(8, 0))
        frame.columnconfigure(5, weight=1)
        frame.rowconfigure(14, weight=1)

    def _path_row(self, frame, row, label, variable, command, save=False):
        ttk.Label(frame, text=label).grid(row=row, column=0, sticky="w", pady=3)
        ttk.Entry(frame, textvariable=variable, width=76).grid(row=row, column=1, columnspan=2, sticky="ew", pady=3)
        ttk.Button(frame, text="Browse", command=command).grid(row=row, column=3, sticky="e", padx=(6, 0))

    def choose_root(self):
        selected = filedialog.askdirectory(title="Choose locally synced Dropbox video folder")
        if selected:
            root = Path(selected)
            self.root_var.set(str(root))
            self.output_var.set(str(root / f"ventilation_results_{datetime.now():%Y%m%d_%H%M%S}.csv"))

    def choose_output(self):
        selected = filedialog.asksaveasfilename(defaultextension=".csv", filetypes=[("CSV", "*.csv")])
        if selected:
            self.output_var.set(selected)

    def choose_roi(self, species: str):
        root = Path(self.root_var.get())
        videos = [
            path for path in root.rglob("*")
            if path.suffix.lower() in VIDEO_EXTENSIONS
            and safe_metadata(path.relative_to(root), root.name)["species"] == species
        ] if root.is_dir() else []
        if not videos:
            messagebox.showerror("No matching videos", f"No {species.replace('_', ' ')} video was found. Check its folder or filename label.")
            return
        cap = cv2.VideoCapture(str(videos[0]))
        ok, frame = cap.read()
        cap.release()
        if not ok:
            messagebox.showerror("Cannot open video", f"Could not read {videos[0].name}")
            return
        target = "buccal or gill movement" if species == "epaulette_shark" else "spiracle or branchial movement"
        x, y, width, height = cv2.selectROI(f"Draw {target} ROI, then press Enter", frame, fromCenter=False, showCrosshair=True)
        cv2.destroyAllWindows()
        if width <= 0 or height <= 0:
            return
        image_height, image_width = frame.shape[:2]
        values = (y / image_height, (y + height) / image_height, x / image_width, (x + width) / image_width)
        for variable, value in zip(self.roi_vars[species], values):
            variable.set(f"{value:.4f}")

    def settings(self) -> Settings:
        rois = {species: tuple(float(variable.get()) for variable in variables) for species, variables in self.roi_vars.items()}
        for species, roi in rois.items():
            if not (0 <= roi[0] < roi[1] <= 1 and 0 <= roi[2] < roi[3] <= 1):
                raise ValueError(f"{species} ROI values must fall between 0 and 1, with each end greater than its start")
        minimum, maximum = float(self.min_bpm_var.get()), float(self.max_bpm_var.get())
        if minimum <= 0 or maximum <= minimum:
            raise ValueError("Ventilation range must have positive minimum and larger maximum")
        return Settings(rois=rois, min_bpm=minimum, max_bpm=maximum, prominence=float(self.prominence_var.get()), diagnostics=self.diagnostics_var.get())

    def start(self):
        try:
            root = Path(self.root_var.get()).expanduser()
            output = Path(self.output_var.get()).expanduser()
            settings = self.settings()
            if not root.is_dir():
                raise ValueError("Choose an existing video root folder")
            videos = sorted(p for p in root.rglob("*") if p.is_file() and p.suffix.lower() in VIDEO_EXTENSIONS)
            if not videos:
                raise ValueError("No supported video files were found")
        except Exception as exc:
            messagebox.showerror("Cannot start", str(exc))
            return
        self.log.delete("1.0", "end")
        self.progress.configure(value=0, maximum=len(videos))
        threading.Thread(target=self.worker, args=(root, output, settings, videos), daemon=True).start()

    def worker(self, root: Path, output: Path, settings: Settings, videos: list[Path]):
        diagnostics = output.parent / "ventilation_diagnostics" if settings.diagnostics else None
        rows = []
        for index, video in enumerate(videos, 1):
            self.queue.put(("log", f"{index}/{len(videos)}  {video.relative_to(root)}\n"))
            row = analyse_video(video, root, settings, diagnostics)
            rows.append(row)
            self.queue.put(("progress", index))
            if row["status"] == "failed":
                self.queue.put(("log", f"  Failed: {row['error']}\n"))
            else:
                self.queue.put(("log", f"  {row['ventilation_bpm']} breaths per minute, {row['method']}\n"))
        fields = sorted({key for row in rows for key in row})
        output.parent.mkdir(parents=True, exist_ok=True)
        with output.open("w", newline="", encoding="utf8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)
        self.queue.put(("done", f"Finished {len(rows)} videos. Results: {output}"))

    def _drain_queue(self):
        try:
            while True:
                event, value = self.queue.get_nowait()
                if event == "log":
                    self.log.insert("end", value)
                    self.log.see("end")
                elif event == "progress":
                    self.progress.configure(value=value)
                elif event == "done":
                    self.log.insert("end", "\n" + value + "\n")
                    self.log.see("end")
                    messagebox.showinfo("Batch analysis complete", value)
        except queue.Empty:
            pass
        self.after(150, self._drain_queue)


if __name__ == "__main__":
    VentilationGUI().mainloop()
