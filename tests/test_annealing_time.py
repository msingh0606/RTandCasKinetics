import unittest
import numpy as np
from primer_func import calculate_annealing_time  # Ensure this module is correctly imported

class TestCalculateAnnealingTime(unittest.TestCase):

    def test_basic_annealing_time(self):
        """Test with standard input parameters."""
        P0 = 50e-9  # Primer concentration (50 nM)
        D = 5e-9    # DNA concentration (5 nM)
        k = 1e6     # Reaction rate constant
        target_percentage = 95

        t_anneal, time, P = calculate_annealing_time(P0, D, k, target_percentage)
        
        # Expected time to anneal based on manual calculation
        expected_t_anneal = -np.log((1 - target_percentage / 100)) / (k * D)
        
        self.assertAlmostEqual(t_anneal, expected_t_anneal, places=6)
        self.assertTrue(np.all(P <= P0))  # Ensure primer concentration never exceeds initial value
        self.assertAlmostEqual(P[-1], P0 * np.exp(-k * D * time[-1]), places=6)  # Verify the final value matches equation

    def test_edge_case_full_annealing(self):
        """Test for 100% annealing."""
        P0 = 50e-9
        D = 5e-9
        k = 1e6
        target_percentage = 100

        with self.assertRaises(ValueError):  # Full annealing (100%) leads to P_target = 0, invalid log
            calculate_annealing_time(P0, D, k, target_percentage)

    def test_edge_case_no_annealing(self):
        """Test for 0% annealing."""
        P0 = 50e-9
        D = 5e-9
        k = 1e6
        target_percentage = 0

        t_anneal, time, P = calculate_annealing_time(P0, D, k, target_percentage)
        self.assertEqual(t_anneal, 0)  # No time is required for 0% annealing
        self.assertTrue(np.allclose(P, P0))  # Primer concentration should remain constant

    def test_invalid_inputs(self):
        """Test invalid inputs like negative concentrations."""
        with self.assertRaises(ValueError):
            calculate_annealing_time(-50e-9, 5e-9, 1e6, 95)  # Negative primer concentration

        with self.assertRaises(ValueError):
            calculate_annealing_time(50e-9, -5e-9, 1e6, 95)  # Negative DNA concentration

        with self.assertRaises(ValueError):
            calculate_annealing_time(50e-9, 5e-9, -1e6, 95)  # Negative reaction rate constant

        with self.assertRaises(ValueError):
            calculate_annealing_time(50e-9, 5e-9, 1e6, -10)  # Negative target percentage

        with self.assertRaises(ValueError):
            calculate_annealing_time(50e-9, 5e-9, 1e6, 100)  # Target percentage = 100%

    def test_output_dimensions(self):
        """Ensure output dimensions match expectations."""
        P0 = 50e-9
        D = 5e-9
        k = 1e6
        target_percentage = 95

        t_anneal, time, P = calculate_annealing_time(P0, D, k, target_percentage)
        
        self.assertEqual(len(time), len(P))  # Ensure time and P arrays have matching lengths
        self.assertTrue(len(time) > 0)       # Ensure arrays are non-empty

if __name__ == "__main__":
    unittest.main()
