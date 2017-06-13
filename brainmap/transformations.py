import numpy as np


def translation_M(v: np.ndarray) -> np.ndarray:
    M = np.eye(4)
    M[:-1, -1] = v
    return M


def scale_M(s: np.ndarray) -> np.ndarray:
    M = np.eye(4)
    M[np.diag_indices(3)] = s
    return M


def shear_M(h: np.ndarray) -> np.ndarray:
    M = np.eye(4)
    M[np.repeat(np.arange(3), 2), np.tile(np.arange(3), 2)] = h
    return M


def rotation_x_M(thetax: float) -> np.ndarray:
    M = np.eye(4)
    M[1, 1] = np.cos(thetax)
    M[1, 2] = -np.sin(thetax)
    M[2, 1] = np.sin(thetax)
    M[2, 2] = np.cos(thetax)
    return M


def rotation_y_M(thetay: float) -> np.ndarray:
    M = np.eye(4)
    M[0, 0] = np.cos(thetay)
    M[0, 2] = np.sin(thetay)
    M[2, 0] = -np.sin(thetay)
    M[2, 2] = np.cos(thetay)
    return M


def rotation_z_M(thetaz: float) -> np.ndarray:
    M = np.eye(4)
    M[0, 0] = np.cos(thetaz)
    M[0, 1] = -np.sin(thetaz)
    M[1, 0] = np.sin(thetaz)
    M[1, 1] = np.cos(thetaz)
    return M


def omogeneous_coordinates(x: np.ndarray) -> np.ndarray:
    if x.shape[0] != 3:
        tmp = np.transpose(x)
    else:
        tmp = x
    return np.vstack([tmp, np.ones(tmp.shape[1])])


def apply_transform(x: np.ndarray, M: np.ndarray) -> np.ndarray:
    tmp = omogeneous_coordinates(x)
    tmp = M.dot(tmp)
    return tmp[:-1, :].T
