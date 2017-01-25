
load test_tourbillon_matlab
U_lambda_matlab = U_lambda;
U_energy_matlab = U_mode_energy;
U_frequency_matlab = U_mode_frequency;
load test_tourbillon_python
U_lambda_python = U_lambda;
U_energy_python = U_mode_energy;
U_frequency_python = U_mode_frequency;
for i=1:size(U_lambda_matlab,1)
    erreur_energy(i) = U_energy_matlab(i)- U_energy_python(i);
    erreur_frequency(i) = U_frequency_matlab(i) - U_frequency_python(i);
end
mean_frequence = mean(erreur_frequency)
mean_energy = mean(erreur_energy)
var_frequence = var(erreur_frequency)
var_energy = var(erreur_energy)