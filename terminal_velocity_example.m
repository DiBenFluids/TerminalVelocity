% Examples


%% positively buoyant disc

[w_sph , w_shape ] = terminal_velocity_master('majoraxis',.01,'minoraxis',.0008,'rho_p',950,'spheroid','oblate')


%% negatively buoyant sphere

[w_sph , w_shape ] = terminal_velocity_master('diameter',.007,'rho_p',1010)

[w_sph , w_shape ] = terminal_velocity_master('diameter',.007,'rho_p',1010 , 'method', 'stokes')

[w_sph , w_shape ] = terminal_velocity_master('diameter',.007,'rho_p',1010 , 'aspectratio', 10)



