import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as st
import networkx as nx
import scipy
from scipy import signal
import cmath

def construct_graph(p,  n_conn_nodes):
    g =  nx.watts_strogatz_graph(n_conn_nodes, 3, 0.5)  #param d modificare
    # adj = (nx.adjacency_matrix(g)).toarray()

    sel_sensors = sample_sensors(p, n_conn_nodes)
    d = dict(zip([i for i in g.nodes()], [int(sel_sensors[j]) for j in range(len(sel_sensors))]))
    g = nx.relabel_nodes(g, d)
    return g

def sample_sensors(p, n_conn_nodes):
    n = np.random.randint(0, p, size = n_conn_nodes)  #65 da modificare
    return n


def construct_theta_subblocks(m):
    block = np.eye(m)*0.2
    return block


def construct_precision_theta(p,m, g):
    delta = 1
    theta = np.eye((m*p))* delta

    for i in g.edges:
        # print(i)
        theta[i[0]*m:(i[0]+1)*m, i[1]*m:(i[1]+1)*m] = construct_theta_subblocks(m)
        theta[i[1]*m:(i[1]+1)*m, i[0]*m:(i[0]+1)*m] = construct_theta_subblocks(m) 

    while delta <= np.max(np.sum(abs(theta -np.eye(m*p)), axis=0)):
        # print(delta, np.max(np.sum(abs(theta -np.eye(m*p)), axis=0)))
        delta += 1

    theta =  theta + np.eye(m*p)*(delta-1)
    print(delta,np.max(np.sum(abs(theta -np.eye(m*p)), axis=0)))
    
    return theta


def define_func_score(p,m, mu, sigma):
    # mu \in R^{pm} ?
    # sigma = theta^(-1), theta \in R^{pm x om}   
    csi = np.random.multivariate_normal(np.zeros(m*p), sigma).reshape((m,p))

    #NB: di può decidere di richiamare la funzione precision matrix qui per creare theta ed invertirla (bisogna vedere se è utile avere theta nel codice, forse sì)
    return csi


def create_fourier_basis(t, freq_band):
    t = np.array(t)
    fb = np.ones((len(freq_band)*2,len(t)))*(-3)
    k=0
    for f in freq_band:
        fb[k, :] = np.sqrt(f)*np.cos(2*np.pi*f*t)
        fb[k+1, :] = np.sqrt(f)*np.sin(2*np.pi*f*t)
        k = k +2
    return fb    


def func_data(t, freq_band, p, mu, sigma):
    fb_vec = create_fourier_basis(t,freq_band).T   #NB non ottimale la generazione della base a ogni dato funzionale generato 
    # print(fb_vec.shape)
    csi_vec = define_func_score(p,len(freq_band)*2,mu,sigma)
    # print(csi_vec.shape)
    sm = np.matmul(fb_vec, (csi_vec))
    return sm.T, csi_vec


def fdata_sim(n1, n2, p, m,  freq_band, time, sigma_pop, sigma_group1):
    data = np.zeros((p, len(time), (n1+n2)))
    csi_vec_pop = np.zeros((m, p, (n1+n2)))
    csi_vec_group = np.zeros((m, p, (n1)))
    for i in range(n1+n2):
        data[:,:, i],  csi_vec_pop[:,:,i] = func_data(time, freq_band, p,  0, sigma_pop)  
        if i < n1:
            fd = func_data(time, freq_band, p,  0, sigma_group1)
            data[:, :, i] += fd[0]
            csi_vec_group[:,:,i]  = fd[1]

    return data, csi_vec_pop, csi_vec_group 


######################################################################
p = 64; 
freq = [9,10,11,12];
m= len(freq)*2; 
n_nodes_active_pop = 15; 
n_nodes_diff_group = 7
T =100;  nt= 10000; deltaT = 1/nt; 
n_sub1 = 24; n_sub2 = 26; n_sub = n_sub1 + n_sub2
time = np.linspace(0, T, nt)

par = {}
par['p'] = p
par['m'] = m
par['freq'] = freq
par['n_nodes_active_pop'] = n_nodes_active_pop
par['n_nodes_diff_group'] = n_nodes_diff_group
par['T'] = T
par['nt'] = nt
par['deltaT'] = deltaT
par['n_sub1'] = n_sub1
par['n_sub2'] = n_sub2
par['n_sub'] = n_sub

np.save("conditional_neurofgm/sim_data_50sub/parameters", par)


g = construct_graph(p, n_nodes_active_pop)
np.savetxt('conditional_neurofgm/sim_data_50sub/adj_pop.txt', np.array(list(g.edges)).astype(int), fmt='%i')

# print(a)
nx.draw_networkx(g, with_labels = True)
plt.savefig('conditional_neurofgm/sim_data_50sub/graph_pop.png')
# plt.show()
plt.close()



tp = construct_precision_theta(p,m,g)
np.savetxt('conditional_neurofgm/sim_data_50sub/theta_pop.txt', tp)
plt.imshow(tp)
plt.savefig('conditional_neurofgm/sim_data_50sub/theta_pop.png')
# plt.show()
plt.close()


# while np.all(np.linalg.eigvals(tp) > 0) != True:
print(np.all(np.linalg.eigvals(tp) > 0))
sigma = np.linalg.inv(tp)

# plt.imshow(sigma)
# plt.show()
# plt.close()


g_diff = construct_graph(p, n_nodes_diff_group)
while len(np.intersect1d(list(g.nodes()), list(g_diff.nodes()))) != 0:
    g_diff = construct_graph(p, n_nodes_diff_group)

np.savetxt('conditional_neurofgm/sim_data_50sub/adj_group.txt', np.array(list(g_diff.edges)).astype(int), fmt='%i')

nx.draw_networkx(g_diff, with_labels = True, edge_color = 'r')
plt.savefig('conditional_neurofgm/sim_data_50sub/graph_group.png')
# plt.show()
plt.close()
   

tp_g1 = construct_precision_theta(p,m,g_diff)
np.savetxt('conditional_neurofgm/sim_data_50sub/theta_group.txt', tp_g1)
plt.imshow(tp_g1)
plt.savefig('conditional_neurofgm/sim_data_50sub/theta_group.png')
# plt.show()
plt.close()

print(np.all(np.linalg.eigvals(tp_g1) > 0))
sigma_g1 = np.linalg.inv(tp_g1)
# plt.imshow(sigma_g1)
# plt.show()
# plt.close()


d, csi_pop, csi_group = fdata_sim(n_sub1,n_sub2,p,m , freq,time,sigma, sigma_g1)   #da risistemare

# d, csi = func_data(time, freq, p, 0, sigma)
print('csi pop', csi_pop.shape)
print('csi group', csi_group.shape)
print('data', d.shape)

for i in range(n_sub):
    np.save("conditional_neurofgm/sim_data_50sub/scores_p_"+str(i), csi_pop[:,:,i])
    np.save("conditional_neurofgm/sim_data_50sub/fd_data_"+str(i), d[:,:, i])
    if i < n_sub1:
            np.save("conditional_neurofgm/sim_data_50sub/scores_g_"+str(i), csi_group[:,:,i])