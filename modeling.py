from __future__ import annotations

import math
import statistics
from dataclasses import dataclass


def solve_linear(matrix: list[list[float]], values: list[float]) -> list[float]:
    size = len(values)
    a = [matrix[i][:] + [values[i]] for i in range(size)]
    for col in range(size):
        pivot = max(range(col, size), key=lambda row: abs(a[row][col]))
        a[col], a[pivot] = a[pivot], a[col]
        divisor = a[col][col]
        if abs(divisor) < 1e-10:
            raise ValueError("模型矩阵不可求解")
        for j in range(col, size + 1):
            a[col][j] /= divisor
        for row in range(size):
            if row == col:
                continue
            factor = a[row][col]
            for j in range(col, size + 1):
                a[row][j] -= factor * a[col][j]
    return [a[i][size] for i in range(size)]


@dataclass
class CsvRidgeModel:
    feature_names: list[str]
    medians: list[float]
    scales: list[float]
    coefficients: list[float]

    def vector(self, row: dict) -> list[float]:
        result = [1.0]
        for name, median, scale in zip(self.feature_names, self.medians, self.scales):
            value = row.get(name)
            result.extend([((median if value is None else value) - median) / scale, float(value is None)])
        return result

    def predict(self, row: dict) -> float:
        return max(0.0, sum(value * coefficient for value, coefficient in zip(self.vector(row), self.coefficients)))


def fit_csv_ridge(samples: list[dict], target: str, feature_names: list[str]) -> CsvRidgeModel:
    labeled = [row for row in samples if row.get(target) is not None]
    if len(labeled) < 8:
        raise ValueError(f"{target} 的训练标签不足 8 个")
    names, medians, scales = [], [], []
    for name in feature_names:
        observed = [row[name] for row in labeled if row.get(name) is not None]
        if not observed:
            continue
        median = statistics.median(observed)
        scale = max(1e-6, math.sqrt(sum((value - median) ** 2 for value in observed) / len(observed)))
        names.append(name)
        medians.append(median)
        scales.append(scale)
    model = CsvRidgeModel(names, medians, scales, [])
    size = 1 + 2 * len(names)
    gram = [[0.0] * size for _ in range(size)]
    cross = [0.0] * size
    for row in labeled:
        vector = model.vector(row)
        for i in range(size):
            cross[i] += vector[i] * row[target]
            for j in range(size):
                gram[i][j] += vector[i] * vector[j]
    for i in range(1, size):
        gram[i][i] += 1.0
    model.coefficients = solve_linear(gram, cross)
    return model
