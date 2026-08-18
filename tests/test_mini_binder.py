import unittest

import torch

from affinity_benchmark.educational.mini_binder import (
    MiniEquivariantDenoiser,
    MiniProteinMPNN,
    UNKNOWN_TOKEN,
    apply_rotation,
    cosine_schedule,
    interface_score,
    make_synthetic_complex,
    q_sample,
    random_rotation,
    sample_with_score_guidance,
)


class MiniBinderTeachingModelTests(unittest.TestCase):
    def setUp(self):
        torch.manual_seed(7)

    def test_schedule_and_forward_diffusion_shapes(self):
        data = make_synthetic_complex(batch_size=3, rotate=False)
        schedule = cosine_schedule(20)
        timestep = torch.tensor([0, 10, 19])
        noisy, noise = q_sample(data.binder, timestep, schedule["alpha_bar"])
        self.assertEqual(noisy.shape, data.binder.shape)
        self.assertEqual(noise.shape, data.binder.shape)
        self.assertTrue(torch.all(schedule["alpha_bar"][1:] < schedule["alpha_bar"][:-1]))

    def test_denoiser_is_rotation_and_translation_equivariant(self):
        data = make_synthetic_complex(batch_size=2, rotate=False)
        model = MiniEquivariantDenoiser(hidden_dim=24, layers=2, max_steps=20)
        timestep = torch.tensor([4, 11])
        noisy = data.binder + 0.3 * torch.randn_like(data.binder)
        reference = model(noisy, data.target, timestep)

        rotation = random_rotation(2)
        translation = torch.tensor([[[2.0, -3.0, 1.5]], [[-1.0, 4.0, 0.25]]])
        rotated_noisy = apply_rotation(noisy, rotation) + translation
        rotated_target = apply_rotation(data.target, rotation) + translation
        transformed = model(rotated_noisy, rotated_target, timestep)
        expected = apply_rotation(reference, rotation) + translation
        self.assertTrue(torch.allclose(transformed, expected, atol=2e-5, rtol=2e-5))

    def test_sequence_model_output_shape(self):
        data = make_synthetic_complex(batch_size=4)
        model = MiniProteinMPNN(hidden_dim=24, layers=2)
        partial = torch.full_like(data.binder_tokens, UNKNOWN_TOKEN)
        logits = model(data.binder, data.target, data.target_tokens, partial)
        self.assertEqual(logits.shape, (4, data.binder.shape[1], 6))

    def test_score_guided_decoder_returns_complete_sequences(self):
        data = make_synthetic_complex(batch_size=2)
        model = MiniProteinMPNN(hidden_dim=16, layers=1)
        sequence, log_probability = sample_with_score_guidance(
            model,
            data.binder,
            data.target,
            data.target_tokens,
            beta=0.5,
        )
        self.assertEqual(sequence.shape, data.binder_tokens.shape)
        self.assertEqual(log_probability.shape, (2,))
        self.assertFalse(torch.any(sequence == UNKNOWN_TOKEN))

    def test_complementary_sequence_scores_better_than_charged_mismatch(self):
        data = make_synthetic_complex(batch_size=1, rotate=False)
        complementary = interface_score(
            data.binder, data.target, data.binder_tokens, data.target_tokens
        )
        mismatched = torch.full_like(data.binder_tokens, 3)
        poor = interface_score(data.binder, data.target, mismatched, data.target_tokens)
        self.assertGreater(complementary.item(), poor.item())


if __name__ == "__main__":
    unittest.main()
