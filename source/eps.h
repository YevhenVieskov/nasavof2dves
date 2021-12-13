      eps_E = 0
      eps_W = 0
      eps_N = 0
      eps_S = 0
      IF (.NOT.(FLAG(i+1,j) < C_F)) eps_E = 1
      IF (.NOT.(FLAG(i-1,j) < C_F)) eps_W = 1
      IF (.NOT.(FLAG(i,j+1) < C_F)) eps_N = 1
      IF (.NOT.(FLAG(i,j-1) < C_F)) eps_S = 1
