import unittest

import numpy as np

from affinity_benchmark.metrics.sequence_ensemble import (
    pairwise_identity_matrix,
    position_frequencies,
    sequence_array,
    shannon_entropy,
)


class SequenceEnsembleMetricTests(unittest.TestCase):
    def test_pairwise_identity_uses_aligned_positions(self):
        observed = pairwise_identity_matrix(["AAA", "ACA", "ACC"])
        expected = np.array(
            [
                [1.0, 2.0 / 3.0, 1.0 / 3.0],
                [2.0 / 3.0, 1.0, 2.0 / 3.0],
                [1.0 / 3.0, 2.0 / 3.0, 1.0],
            ]
        )
        np.testing.assert_allclose(observed, expected)

    def test_position_frequencies_and_entropy(self):
        frequencies = position_frequencies(["AAA", "ACA", "ACC"], alphabet="AC")
        np.testing.assert_allclose(
            frequencies,
            np.array([[1.0, 0.0], [1.0 / 3.0, 2.0 / 3.0], [2.0 / 3.0, 1.0 / 3.0]]),
        )
        entropy = shannon_entropy(frequencies)
        expected_mixed = -(1.0 / 3.0) * np.log(1.0 / 3.0) - (2.0 / 3.0) * np.log(
            2.0 / 3.0
        )
        np.testing.assert_allclose(entropy, [0.0, expected_mixed, expected_mixed])

    def test_rejects_unequal_lengths_and_unknown_residues(self):
        with self.assertRaisesRegex(ValueError, "common length"):
            sequence_array(["AA", "AAA"])
        with self.assertRaisesRegex(ValueError, "outside the alphabet"):
            sequence_array(["AAZ"])

    def test_rejects_non_normalized_frequencies(self):
        with self.assertRaisesRegex(ValueError, "sum to one"):
            shannon_entropy(np.array([[0.2, 0.2]]))


if __name__ == "__main__":
    unittest.main()
