import mujoco
from pathlib import Path
MODEL_PATH = Path(__file__).with_name("planar_leg.xml")
model = mujoco.MjModel.from_xml_path(str(MODEL_PATH))
data = mujoco.MjData(model)

print("Has data.qM:", hasattr(data, "qM"))
print("Has data.qfrc_bias:", hasattr(data, "qfrc_bias"))
print("Type:", type(data))