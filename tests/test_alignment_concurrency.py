import time

from src.alignment import iter_bounded_results


def test_bounded_results_limits_pending_work_and_runs_each_task_once():
    active = 0
    peak = 0
    seen = []

    def worker(task):
        nonlocal active, peak
        active += 1
        peak = max(peak, active)
        time.sleep(0.002)
        seen.append(task)
        active -= 1
        return task * 2

    results = list(iter_bounded_results(range(12), worker, workers=2, max_pending=3))
    assert sorted(value for _, value in results) == [index * 2 for index in range(12)]
    assert sorted(seen) == list(range(12))
    assert peak <= 2
