"""design.h5 / result.h5 的读写契约。"""

import h5py
import numpy as np
import pytest

import dataio


# ------------------------------------------------------------------ design

def test_design_roundtrip(tmp_path):
    path = str(tmp_path / "design.h5")
    in_P = np.arange(24, dtype=float).reshape(12, 2)
    ex_P = in_P + 1.0
    ports = dict(
        port_id=np.array([[3, 4, 5, 6]]),
        port_is_analytic=np.array([[1, 1, 1, 1]]),
        port_coord=np.zeros((3, 4)),
        port_x_axis=np.ones((3, 4)),
        incident_on=np.array([[1, 0, 0, 0]]),
        port_other_params=np.full((8, 4), 0.5),
    )
    dataio.write_design(path, in_P=in_P, in_V=in_P, ex_P=ex_P, ex_V=ex_P, **ports)

    back = dataio.read_design(path)
    assert np.array_equal(back["in_Bezier_P"], in_P)
    assert np.array_equal(back["ex_Bezier_P"], ex_P)
    assert np.array_equal(back["ports"]["port_id"], ports["port_id"])
    assert np.array_equal(back["ports"]["port_other_params"], ports["port_other_params"])
    assert back["meta"]["format_version"] == dataio.DESIGN_FORMAT_VERSION
    # 整数表应以整型存，实数表以 float64 存
    assert back["ports"]["port_id"].dtype.kind == "i"
    # 没给 node_t / dofs 时不应凭空捏造
    assert back["node_t"] is None
    assert back["dofs"] is None


def test_design_node_t_roundtrip(tmp_path):
    """node_t 用 NaN 表示"不在优化边界上"，必须原样往返。"""
    path = str(tmp_path / "design.h5")
    node_t = np.array([np.nan, 0.0, 3.5, np.nan, 11.0])
    dataio.write_design(
        path,
        in_P=np.zeros((3, 2)), in_V=np.zeros((3, 2)),
        ex_P=np.zeros((3, 2)), ex_V=np.zeros((3, 2)),
        port_id=np.zeros((1, 4)), port_is_analytic=np.zeros((1, 4)),
        port_coord=np.zeros((3, 4)), port_x_axis=np.zeros((3, 4)),
        incident_on=np.zeros((1, 4)), port_other_params=np.zeros((8, 4)),
        node_t=node_t,
    )
    back = dataio.read_design(path)["node_t"]
    assert back.shape == node_t.shape
    assert np.array_equal(np.isnan(back), np.isnan(node_t))
    assert np.allclose(back[np.isfinite(back)], node_t[np.isfinite(node_t)])


# ------------------------------------------------------------------ result

def _write_result(path, S, sens_ex, sens_in, ez=None, version=dataio.FORMAT_VERSION):
    with h5py.File(path, "w") as f:
        f.create_group("meta").attrs["format_version"] = version
        f["meta"].attrs["n_process"] = 3
        sc = f.create_group("scalars")
        sc.create_dataset("S_params", data=np.stack([S.real, S.imag], axis=1))
        sc.create_dataset("ex_opt_sens", data=sens_ex)
        sc.create_dataset("in_opt_sens", data=sens_in)
        if ez is not None:
            f.create_group("fields").create_dataset(
                "Ez", data=np.stack([ez.real, ez.imag], axis=1))


def test_read_result(tmp_path):
    path = str(tmp_path / "result.h5")
    S = np.array([1 + 2j, 3 - 4j])
    sens = np.arange(6, dtype=float).reshape(3, 2)
    ez = np.array([1 + 1j, 2 - 2j, 0.5 + 0j])
    _write_result(path, S, sens, -sens, ez=ez)

    res = dataio.read_result(path, want_field=True)
    assert np.allclose(res["S_params"], S)
    assert np.array_equal(res["sens_ex"], sens)
    assert np.array_equal(res["sens_in"], -sens)
    assert np.allclose(res["Ez"], ez)
    assert res["n_process"] == 3

    # 不要求场时不应去碰 /fields
    assert "Ez" not in dataio.read_result(path)


def test_read_result_without_field(tmp_path):
    path = str(tmp_path / "result.h5")
    _write_result(path, np.array([1 + 0j]), np.zeros((2, 2)), np.zeros((2, 2)))
    assert "Ez" not in dataio.read_result(path, want_field=True)


def test_read_result_rejects_old_layout(tmp_path):
    """旧文件（format_version=1，每 rank 一份）要给一句能看懂的话。"""
    path = str(tmp_path / "result.h5")
    with h5py.File(path, "w") as f:
        f.create_group("meta").attrs["format_version"] = 1
        sc = f.create_group("scalars")
        sc.create_dataset("S_params_0", data=np.zeros((4, 2)))
        sc.create_dataset("ex_opt_sens", data=np.zeros((2, 2)))
        sc.create_dataset("in_opt_sens", data=np.zeros((2, 2)))

    with pytest.raises(ValueError, match="格式版本"):
        dataio.read_result(path)


def test_read_result_missing_file(tmp_path):
    with pytest.raises(FileNotFoundError):
        dataio.read_result(str(tmp_path / "nope.h5"))
