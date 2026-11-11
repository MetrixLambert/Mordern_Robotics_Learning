from pathlib import Path
import time
import mujoco
import mujoco.viewer
import numpy as np

xml = Path(__file__).with_name('planar_leg.xml')
model = mujoco.MjModel.from_xml_path(str(xml))
data = mujoco.MjData(model)
mujoco.mj_resetDataKeyframe(model, data, 0)
mujoco.mj_forward(model, data)
print('nq, nv, nu =', model.nq, model.nv, model.nu)
print('joint order:', [mujoco.mj_id2name(model, mujoco.mjtObj.mjOBJ_JOINT, i) for i in range(model.njnt)])
print('initial qpos:', data.qpos.copy())
print('initial foot center:', data.site_xpos[model.site('foot_center').id].copy())

with mujoco.viewer.launch_passive(model, data) as viewer:
    for step in range(3000):
        t0 = time.perf_counter()
        data.ctrl[:] = 0.0  # no control yet: observe gravity and ground contact
        mujoco.mj_step(model, data)
        if step % 250 == 0:
            print(f't={data.time:.2f}s  q={np.round(data.qpos, 4)}  contacts={data.ncon}')
        viewer.sync()
        time.sleep(max(0.0, model.opt.timestep - (time.perf_counter() - t0)))
