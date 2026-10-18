"""C++ 求解器控制：写 ``infile.pot``、调用 ``mpiexec``、校验、读结果。

这一层只负责"把 Case 交给求解器，把 Result 拿回来"，不含任何几何或优化逻辑。
"""

import os

import dataio
import paths

# 默认求解器选项（直接法 LU）
DEFAULT_SOLVER_ARGS = "-pc_type lu"


def infile_lines(case, save_field=False):
    """生成 ``infile.pot`` 的内容。

    键名必须与 C++ 侧 ``WaveguideParams`` / ``main.cpp`` 读的一致。
    ``id_*`` 是 gmsh 里的数字 physical id；给出时 C++ 就不必再按
    physical name 查（名字在格式转换时容易丢）。
    """
    ph = case.physics
    lines = [
        ("mesh_file", case.mesh_file),
        ("out_mesh_file", case.out_mesh_file),
        ("port_id_varname", "port_id"),
        ("port_analytic_varname", "port_is_analytic"),
        ("port_coords_varname", "port_coord"),
        ("port_other_params_varname", "port_other_params"),
        ("port_x_axis_varname", "port_x_axis"),
        ("port_incidence_on_varname", "incident_on"),
        ("start_point_mode", "0"),
        ("end_point_mode", "0"),
        ("design_file", case.design_file),
        ("result_file", case.result_file),
        ("epsilon_wg", "%f" % ph.epsilon_wg),
        ("epsilon_amb", "%f" % ph.epsilon_amb),
        ("mu_wg", "%f" % ph.mu_wg),
        ("mu_amb", "%f" % ph.mu_amb),
        ("freq", "%f" % ph.freq),
    ]
    ids = case.ids or {}
    for key, name in (("wg", "id_wg"), ("ex_opt", "id_ex_opt"), ("in_opt", "id_in_opt")):
        if key in ids:
            lines.append((name, str(ids[key])))
    if save_field:
        lines.append(("save_field", "1"))
    return ["%s = %s" % pair for pair in lines]


def write_infile(case, save_field=False, filename="infile.pot"):
    with open(filename, "w") as fp:
        for line in infile_lines(case, save_field=save_field):
            print(line, file=fp)
    return filename


def run(case, extra_args=DEFAULT_SOLVER_ARGS, save_field=False, check=True):
    """写配置 → 跑一次正问题 → 返回解析好的结果字典。"""
    write_infile(case, save_field=save_field)

    # 先删掉旧结果，否则求解器中途失败时会读到上一次的数字而毫无提示
    try:
        os.remove(case.result_file)
    except FileNotFoundError:
        pass

    status = paths.run_solver(extra_args)
    if check and status != 0:
        raise RuntimeError(
            "求解器失败（退出码 %s）。请检查 data/ 下的输出。" % status)
    return dataio.read_result(case.result_file, want_field=bool(save_field))


def read_result(case, want_field=False):
    return dataio.read_result(case.result_file, want_field=want_field)
