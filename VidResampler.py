import cv2
from background_utilities import grab_vids

FOLDER = 'E:/GE/GE Ethanol Quadrants' 


def reSampleVideo(input_path, output_path, frame_interval=1000):
    cap = cv2.VideoCapture(input_path)
    if not cap.isOpened():
        print(f"Error: Could not open video file {input_path}")
        return
    width = int(cap.get(cv2.CAP_PROP_FRAME_WIDTH))
    height = int(cap.get(cv2.CAP_PROP_FRAME_HEIGHT))
    fps = cap.get(cv2.CAP_PROP_FPS)
    total_frames = int(cap.get(cv2.CAP_PROP_FRAME_COUNT))
    fourcc = cv2.VideoWriter_fourcc(*'XVID')
    out = cv2.VideoWriter(output_path, fourcc, fps, (width, height))
    print(f"Total frames in source: {total_frames}")
    print(f"Extracting 1 out of every {frame_interval} frames...")
    
    current_frame = 0
    saved_count = 0
    while current_frame < total_frames:
        cap.set(cv2.CAP_PROP_POS_FRAMES, current_frame)
        ret, frame = cap.read()
        if not ret:
            break
        out.write(frame)
        saved_count += 1
        current_frame += frame_interval
    cap.release()
    out.release()
    print(f"Process complete! Saved {saved_count} frames to {output_path}")

def batch_track(folder):
    vids = grab_vids(folder)
    for i in vids:
        OUTPUT_FILE = i + 'output_sampled.avi'
        reSampleVideo(i, OUTPUT_FILE, frame_interval=75)


if __name__ == "__main__":
    batch_track(folder=FOLDER)