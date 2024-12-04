import unittest

1  # Define the function
def calculate_dNTP_probability(base, dNTP_Conc, NRTI_Conc, Kaff):
    if base == "T":
        fraction = dNTP_Conc / (dNTP_Conc + Kaff * NRTI_Conc)
    else:
        fraction = 1.0  # No competition for non-T bases
    return fraction

# Define the test case
class TestCalculateDNTPProbability(unittest.TestCase):
    
    def test_non_T_base(self):
        """Test calculation when base is not 'T'."""
        result = calculate_dNTP_probability('A', 10, 5, 2)
        self.assertEqual(result, 1.0, "Probability should be 1.0 for non-T bases regardless of inputs.")

    def test_T_base_with_non_zero_values(self):
        """Test calculation for base 'T' with valid NRTI_Conc and Kaff."""
        result = calculate_dNTP_probability('T', 10, 5, 2)
        expected_fraction = 10 / (10 + 2 * 5)
        self.assertAlmostEqual(result, expected_fraction, places=6, msg="Incorrect probability for base 'T'.")
    
    def test_edge_case_n_zero(self):
        """Test calculation when n=0 (probability should always be 1)."""
        result = calculate_dNTP_probability('T', 10, 5, 2, 0)
        self.assertEqual(result, 1.0, "Probability should be 1.0 when n=0.")

    def test_edge_case_zero_dNTP_Conc(self):
        """Test calculation when dNTP_Conc=0 for base 'T'."""
        result = calculate_dNTP_probability('T', 0, 5, 2)
        self.assertEqual(result, 0.0, "Probability should be 0.0 when dNTP_Conc is 0.")
    
    def test_large_n(self):
        """Test calculation with a large n to check behavior."""
        result = calculate_dNTP_probability('T', 10, 5, 2)
        expected_fraction = 10 / (10 + 2 * 5)
        self.assertAlmostEqual(result, expected_fraction, places=6, msg="Incorrect probability for large n.")

if __name__ == "__main__":
    unittest.main()
