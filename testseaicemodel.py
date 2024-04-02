import pytest

from model.rhs import apply_boundary_condition, correct_for_brine_movement

import userinput

# TODO: Add fourier time step case

class TestSeaIceModel:
    def test_userinput(self)-> None:
        assert userinput.constants in ["debug", "real"], "constants must be either debug or real"
        assert userinput.geom in [1, 2], "geom must be either 1 or 2"
        assert userinput.iter_max > 0, "iter_max must be greater than 0"
        assert userinput.Stefan in [True, False], "Stefan must be either True or False"
        assert userinput.Buffo in [True, False], "Buffo must be either True or False"
        assert userinput.liq_rel in ["Normal", "FREZCHEM"], "liq_rel must be either Normal or FREZCHEM"
        assert userinput.dz > 0, "dz must be greater than 0"
        assert userinput.bc_condition in ["Dirichlet", "Neumann"], "bc_condition must be either Dirichlet or Neumann"
        assert userinput.T_tol > 0, "T_tol must be greater than 0"
        assert userinput.S_tol > 0, "S_tol must be greater than 0"
        assert userinput.phi_tol > 0, "phi_tol must be greater than 0"
        assert userinput.T_IC in ["Tm_w", "T(S)", "T271.25", "T_Stefan"], "T_IC must be either Tm_w, T(S), T271.25 or T_Stefan"
        assert userinput.S_IC in ["S_linear", "SX"], "S_IC must be either S_linear or SX where X is a salinity value; ex S33 with salinity being 33"
        assert userinput.P_IC in ["P1", "P_Stefan", "P0"], "P_IC must be either P1, P_Stefan or P0"
        assert userinput.top_temp in ["Stefan", "T_const_250", "T_const_260", "T_const_265", "T_W3"], "top_temp must be either Stefan, T_const_250, T_const_260, T_const_265 or T_W3"
        assert userinput.dt > 0, "dt must be greater than 0 and less than {}".format(userinput.fourier_number_timestep())
        assert userinput.topT_stefan == userinput.T_air_Stefan, "topT_stefan must be equal to T_air_Stefan"