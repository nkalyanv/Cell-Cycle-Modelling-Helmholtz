from tracemalloc import start
import scipy.io
import numpy as np
import sys
import random
import copy
import Utils
import seaborn as sns
import matplotlib.pyplot as plt

medium = 'SCGE'

time_of_protein_production = {'YPD' : 24, 'SCD' : 28, 'SCGE' : 35.5}

trackback_time = {'SCGE' : 13, 'SCD' : 10, 'YPD' : 7}
mrna_amount = {'YPD' : 30, 'SCD' : 21, 'SCGE' : 10}
k_pre = {'SCGE' : (0.000979554896913755, -0.0246254114090489), 
        'SCD' : (0.00117358418816900, -0.0337438108854191),
        'YPD' : (0.00128952510678458, -0.0392751771327524)}
k_post = {'SCGE' : (0.000978964053106707,  0, -0.0107700336183466),
        'SCD' : (0.00102472091756443, 0, -0.0134007070518360),
        'YPD' : (0.000875120797613734, 0, -0.00804804723947798)}

#Growth Rate 1-m_growth -> p1(1), c_growth -> p1(2)
m_growth = {'SCGE' : 0.00222019056145497, 'SCD' : 0.00450131417714158, 'YPD' : 0.00423348155253518}
c_growth = {'SCGE' : 0.104453745081658, 'SCD' : 0.213401235570336, 'YPD' : 0.369664105717204}

#Growth Rate 2
m_growth2 = {'SCGE' : 0.00330368790213606, 'SCD' : 0.0123947445945605, 'YPD' : 0.0142057829609836}  
c_growth2 = {'SCGE' : 0.188960420713999, 'SCD' : -0.166518313652299, 'YPD' : -0.152637732413614}

#Parameters File Paths
parameters_file_path = {'SCGE' : '..\Data\HTBData\HTB2_SCGE_parameters_full.mat', 'SCD' : '..\Data\HTBData\HTB2_SCD_parameters_full.mat', 'YPD' : '..\Data\HTBData\HTB2_YPD_parameters_full.mat'}

def extend(arr):
    for i in range(len(arr)):
        c = arr[i]
        arr[i] = np.concatenate([c, c[:, -1][:, None]], axis = 1)
    return arr

def growth_rate1(vol):
    #p1 in matlab code growthrate_dv_code_new
    m = m_growth[medium]
    c = m_growth[medium]
    #return np.exp(m) + c
    return m*vol + c

def growth_rate2(vol):
    #p11 in above file for single line
    m = m_growth2[medium]
    c = c_growth2[medium]
    return m*vol + c


parameters = scipy.io.loadmat(parameters_file_path[medium])
sys.stdout = open("../Debug/CellCycleDebug.txt", "w")

#Simulation Parameters controlling Number of cells: Total cells = cell_num * cell_bin_num.
# These will be spread equally through the g1_bin_range
cell_num = 40
cell_bin_num = 40
g1_bin_range = [10, 70]

#Cell Parameters
G1_lambda = k_pre[medium]
bud_prob_beta_1 = G1_lambda[0]
bud_prob_M0_1 = -G1_lambda[1]/bud_prob_beta_1

G2_lambda_2d= k_post[medium]

#G1 Simulation Parameters
n_cells = cell_bin_num * cell_num #Total number of cells

birth_size_list = np.linspace(start = g1_bin_range[0], stop = g1_bin_range[1], num = cell_bin_num + 1) #Create list of equi-seperated values-(lowerBound, upperBound, numValues)
birth_size_list = birth_size_list[:-1] + np.diff(birth_size_list) #Shift list up by one position

cell_v = []
for i in birth_size_list:
    cell_v.append(i * np.ones(cell_num))
cell_v = (np.concatenate(cell_v))[:, None]

cell_cycle = np.zeros((n_cells, 1))
cell_bud = np.zeros((n_cells, 1))  #begin with no buds
G2_counter = np.zeros((n_cells, 1))  #cells count down G2
mother_daughter = np.zeros((n_cells, 1)) #begin with all daughters; daughter = 0, mother = 1
mother_G1_counter = np.zeros((n_cells, 1)) #all cells start in G1 as daughters
post_start_G1 = np.zeros((n_cells, 1))
post_start_counter = np.zeros((n_cells, 1))
start_size = np.zeros((n_cells, 1))
cell_list_G1 = set(range(0, n_cells))
bud_prob = np.zeros((n_cells, 1))

cell_volumes = (cell_v).tolist()

#G1 Simulation
size_at_birth = cell_v
size_at_start = np.zeros((n_cells, 1))
G1_duration = np.zeros((n_cells, 1))
G1_growth = np.zeros((n_cells, 1))

current_time = 0

cell_protein = np.zeros((n_cells, 1))
protein_synthesis_counter = np.zeros((n_cells, 1))
protein_synthesis_counter[:] = time_of_protein_production[medium]

while(np.all(size_at_start) == False): #Checks if all values are non-zero. That is every cell has reached start.
    
    current_time += 1
    i = current_time - 1

    bud_prob = bud_prob_beta_1*(cell_v[:, i] - bud_prob_M0_1)  
    bud_prob[bud_prob < 0] = 0
    bud_prob[bud_prob > 1] = 1
    bud_prob = 1 - np.exp(-bud_prob)
    
    bud_len = len(bud_prob)
    bud = np.zeros(bud_len)
    for j in range(bud_len):
        bud[j] = random.choices([0., 1.], weights = [1-bud_prob[j], bud_prob[j]])[0]
    
    #New arrays to be appended
    cell_v, cell_cycle, cell_bud, mother_daughter = extend([cell_v, cell_cycle, cell_bud, mother_daughter])

    for j in copy.copy(cell_list_G1):

        dv = growth_rate1(cell_v[j])       
        cell_v[j]= cell_v[j][i] + dv
        cell_volumes[j].append(cell_v[j][i])
        if(bud[j] == 1): #If Cell passes start
            size_at_start[j] = cell_v[j][i]
            G1_duration[j] = current_time
            G1_growth[j] = cell_v[j][i] - size_at_birth[j]
            cell_list_G1.remove(j)

# G2-Simulation

n_cells = cell_num * cell_bin_num
cell_v = size_at_start
cell_cycle = np.zeros((n_cells,1))
cell_bud = np.zeros((n_cells,1))  
G2_counter = np.zeros((n_cells,1))  
mother_daughter = np.zeros((n_cells,1))
post_start_G1 = np.ones((n_cells,1))
start_size = cell_v

current_time = 0
cell_list_G2 = set(range(0, n_cells))
size_at_div = np.zeros((n_cells,1))
growth_post_start = np.zeros((n_cells,1))
duration_post_start = np.zeros((n_cells,1))
duration_post_start_G1 = np.zeros((n_cells,1))
mother_bud_mass_defect = np.zeros((n_cells,1))



while(len(cell_list_G2) > 0):
    current_time += 1
    i = current_time - 1

    cell_v, cell_bud, cell_cycle, post_start_G1, post_start_counter, G2_counter, cell_protein = extend([cell_v, cell_bud, cell_cycle, post_start_G1, post_start_counter, G2_counter, cell_protein])    

    for k in copy.copy(cell_list_G2):
        dv2 = growth_rate2(cell_v[k][i - 1])
        if(cell_cycle[k][i - 1] == 0): #Cell is in G1(Post-start)
            if(post_start_G1[k][i - 1] == 1): #In post-start G1.
                if(protein_synthesis_counter[k] > 0):
                    cell_protein[k][i] = cell_protein[k][i - 1] + (dv2*mrna_amount[medium]/cell_v[k][i])                  
                    protein_synthesis_counter[k] -= 1
                else:
                    cell_protein[k][i] = cell_protein[k][i - 1]
                if(post_start_counter[k][i] > 0):
                    cell_v[k][i] = cell_v[k][i - 1] + dv2
                    cell_cycle[k][i] = 0
                    cell_bud[k][i] = 0
                    post_start_G1[k][i + 1] = 1
                    post_start_counter[k][i] = post_start_counter[k][i - 1] - 1
                elif(post_start_counter[k][i - 1] == 0):
                    
                    mother_bud_mass_defect[k] = 0
                    cell_v[k][i] = cell_v[k][i] + mother_bud_mass_defect[k]
                    cell_bud[k][i] = cell_bud[k][i - 1] + dv2 - mother_bud_mass_defect[k]
                    cell_cycle[k][i] = 1
                    post_start_G1[k][i] = 0
                    post_start_counter[k][i] = 0
                    #Add new protein amount when shifting to G2
                    protein_extra_time = trackback_time[medium]
                    idx = max(0, len(cell_volumes[k]) - protein_extra_time)
                    for m in range(protein_extra_time):
                        z = min(idx + m, len(cell_volumes[k]) - 1)
                        dv2 = growth_rate2(cell_volumes[k][z])
                        cell_protein[k][i] = cell_protein[k][i] + (dv2*mrna_amount[medium]/cell_volumes[k][z])                  

        elif(cell_cycle[k][i - 1] == 1): #Cell is in G2
            div_prob = G2_lambda_2d[0]*cell_bud[k][i - 1] + G2_lambda_2d[1]*start_size[k] + G2_lambda_2d[2]
            div_prob = 1 - np.exp(-div_prob)
            div_prob = max(0, div_prob)
            div_prob = min(1, div_prob)
            
            G2_counter[k][i - 1] = random.choices([0., 1.], weights = [div_prob, 1 - div_prob])[0]
            
            if(protein_synthesis_counter[k] > 0):
                cell_protein[k][i] = cell_protein[k][i - 1] + (dv2*mrna_amount[medium]/cell_v[k][i])
                protein_synthesis_counter[k] -= 1
            else:
                cell_protein[k][i] = cell_protein[k][i - 1]
            if (G2_counter[k][i - 1] > 0): #growth in G2
                cell_v[k][i] = cell_v[k][i - 1] + mother_bud_mass_defect[k] #bud grows, mother does not
                cell_bud[k][i] = cell_bud[k][i - 1] + dv2 - mother_bud_mass_defect[k]; #grow bud
                cell_cycle[k][i] = 1; #Stay in G2
                post_start_G1[k][i] = 0
                post_start_counter[k][i] = 0

            elif (G2_counter[k][i - 1] == 0):
                size_at_div[k] = cell_v[k][i - 1] + cell_bud[k][i - 1]
                growth_post_start[k] = size_at_div[k] - start_size[k]
                duration_post_start[k] = current_time
                cell_list_G2.remove(k)

Utils.fillForward(cell_protein)

final_protein_content = cell_protein[:, -1]
total_growth = G1_growth + growth_post_start
print("Average protein content:", np.mean(final_protein_content))

np.save(medium + 'protein', final_protein_content)
#Plotting Testing
n_bins = 20
[protein_content,error,fin_size] = Utils.KS_bindata_mean_20140916(size_at_div[:, -1], final_protein_content, n_bins)
plt.errorbar(fin_size, protein_content, error, color = "red")
sns.scatterplot(y = final_protein_content, x = size_at_div[:, -1])
plt.title(medium)
plt.xlabel("Final Cell Size")
plt.ylabel("Protein Content")
plt.savefig('../Plots/ProteinvSize' + medium + '.png')

plt.clf()

plot_parameters = Utils.loadmat(parameters_file_path[medium])
totvolSG2M = plot_parameters['pulsedata']['totvolumebirth']
totvolumeendSG2M = plot_parameters['pulsedata']['totvolumeendSG2M']
sns.scatterplot(totvolSG2M, totvolumeendSG2M, color = 'black')

[volumecytokinesis_mean,volumecytokinesis_error,volumebudemerg_binsmean] = Utils.KS_bindata_mean_20140916(totvolSG2M, totvolumeendSG2M, 20, nargs_out = 3)
plt.errorbar(volumebudemerg_binsmean, volumecytokinesis_mean, volumecytokinesis_error, color = "red")

i = ~np.isnan(volumecytokinesis_mean)
p05 = np.polyfit(volumebudemerg_binsmean[i],volumecytokinesis_mean[i],1)
fit05 = np.polyval(p05, volumebudemerg_binsmean[i])
sns.lineplot(volumebudemerg_binsmean[i],fit05,color = 'green')
[y_size_at_division,error,x_size_at_budemg] = Utils.KS_bindata_mean_20140916(size_at_birth, size_at_div, n_bins, nargs_out = 3)
plt.errorbar(x_size_at_budemg, y_size_at_division, error, color = "blue")
plt.title(medium)
plt.xlabel("volume at birth [fl]")
plt.ylabel("volume at division [fl]")
plt.savefig('..\Plots\model_' + medium + '.png')