"""Run a single blind fragment-anchor trajectory fork from ordinary docking inputs."""
from pathlib import Path
import importlib.util
import sys

from anchors import place

FIELD_ENGINE = Path(__file__).resolve().parents[2] / "cavity_field_rank/engine"
sys.path.insert(0, str(FIELD_ENGINE))
spec = importlib.util.spec_from_file_location("field_live", FIELD_ENGINE / "live.py")
field_live = importlib.util.module_from_spec(spec)
spec.loader.exec_module(field_live)

if __name__ == "__main__":
    field_live.main(generator=place, arms=("whole", "fragment", "sham"))
