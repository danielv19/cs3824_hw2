import math
import random
from time import process_time
import numpy as np
import networkx as nx

#given a key of a node, its cluster is given
Cluster_Dict = {}
#given a key of a cluster, a list of the node within the cluster is given
Cluster_List_Dict =  {}
#given a key of a node, give its neighbors
Neighbors = {}
#num of edges in graph
edges = 0
#output file
outfile = open(f"output_leiden_gammas_k15.txt",'w')
#dict of cluster_degree_sum per cluster
cluster_degree_sum_dict ={}
#dict of edges within cluster
cluster_edges_dict = {}
#dict of edge contribution within cluster
cluster_edge_contribution = {}
#edges within a cluster
global_edge_contr = [0, []]
#modularity when we exit louvain/leiden
exiting_modularity = 0
#len of paritition
max_partition_len = 0
#total change
total_change_mnf = 0

#bypass delta function usage by only looking at nodes within the cluster
def modularity(gamma, partition):
  pas_in = len(cluster_edge_contribution) == 0

  qc_sum = 0
  #if we pass the partition, get_cluster, then lets take the modularity of each cluster
  for cluster in partition:
    qc_sum += cluster_modularity(gamma, partition, cluster)
  
  #initialize the edge contribution dict
  if pas_in:
    for u in Cluster_Dict:
      cluster_edge_contribution[u] = [0,[]]
      cluster_edge_contribution[u][0] = edges_in_cluster(partition,cluster)
  
  return qc_sum

def cluster_modularity(gamma, partition, cluster):
  global cluster_degree_sum_dict
  
  qc_sum = 0
  # to keep track of whether we have seen the pair (u,u). I wish I had a better way of implementing this check.
  seen = []
  pass_in = False
  if (len(cluster_degree_sum_dict) != len (partition)) and len(partition[cluster]) > 0:
    cluster_degree_sum_dict[cluster] = 0
    pass_in = (len(cluster_degree_sum_dict) == len (Cluster_List_Dict))

  
  for u in partition[cluster]:
    if (len(cluster_degree_sum_dict) != len(partition) and len(partition[cluster]) > 0) or (pass_in and len(partition[cluster]) > 0):
      cluster_degree_sum_dict[cluster] += len(Neighbors[u])
      pass_in = False
    
    for v in partition[cluster]:
      if f"{u}{v}" in seen:
        # don't process any edge twice
        next
      else:
        # record that we have seen the pair (u, u)
        seen.append(f'{u}{v}')
      #length of the array containing u's neighbors = degree
      temp = (len(Neighbors[u])) * (len(Neighbors[v])) # d(u) * d(v)
      temp = float(temp)/float(2*edges) # d(u)d(v) / 2m
      #a(u,v) returns 1 if edge (u,v) exists, 0 else
      qc_sum += float(a(u,v) - temp) # a(u,v) - ( d(u)d(v) / 2m )
  qc_sum = qc_sum * (float(1)/float(edges*2))# 1/2m * SUM(seen above)
  return qc_sum * gamma

def change_in_modularity(i,cluster,gamma,get_cluster):
  global global_edge_contr
  global_edge_contr = [0,[]]

  change_sum = 0.0
  for x in [get_cluster[i],cluster]:
    cluster_degree_sum = 0.0
    cluster_degree_sum += cluster_degree_sum_dict[x]

    #act as if i is already in C(j)
    if x == cluster and get_cluster[i] != cluster:
      cluster_degree_sum += (len(Neighbors[i]))
    
    #d(i)/m^2 * Sum of d(u) for u in cluster passed in
    fraction = float(2*edges*edges)
    fraction = float(len(Neighbors[i])) / fraction
    cluster_degree_sum *= (gamma * (fraction)) 
    #1/m * d(i,cluster)
    degree_i = 0.0
    if(x == get_cluster[i]):
      degree_i = float(cluster_edge_contribution[i][0])/float(edges)  
    else: 
      global_edge_contr[0] = cluster_edges(i, x, get_cluster)
      degree_i = float(global_edge_contr[0])/float(edges)
    
    #taking the sum
    if (x != get_cluster[i]):
      change_sum += float(degree_i - cluster_degree_sum)
    else:
      #if we are moving the only node out of a cluster, then the modularity will be 0
      change_sum += float(cluster_degree_sum - degree_i)

  return change_sum

def edges_in_cluster(partition,cluster):
  global cluster_edges_dict
  global cluster_edge_contribution
  edge_sum = 0
  seen = []
  pass_in = False
  if (len(cluster_edges_dict) != len (partition)) and len(partition[cluster]) > 0:
    cluster_edges_dict[cluster] = 0
    pass_in = (len(cluster_edges_dict) == len (partition))
  else:
    return cluster_edges_dict[cluster]
    
  for u in partition[cluster]:
    for v in partition[cluster]:
      if f"{u}{v}" not in seen:
        seen.append(f'{u}{v}')
        seen.append(f'{v}{u}')
        if u in Neighbors[v]:
          if ((len(cluster_edges_dict) != len(partition)) or pass_in) and len(partition[cluster]) > 0:
            cluster_edge_contribution[u][0] += 1
            cluster_edge_contribution[v][0] += 1
            cluster_edge_contribution[v][1].append(u)
            cluster_edge_contribution[u][1].append(v)
          edge_sum += 1
  
  if ((len(cluster_edges_dict) != len(partition)) or pass_in) and len(partition[cluster]) > 0:
    cluster_edges_dict[cluster] += edge_sum
    pass_in = False

  return edge_sum

def change_in_cpm(i,cluster,gamma,partition,get_cluster):
  global global_edge_contr
  change_sum = 0.0
  global_edge_contr = [0,[]]

  for x in [get_cluster[i],cluster]:
    original = cluster_edges_dict[x] - float(gamma*(math.comb(len(partition[x]),2)))
    change_sum -= original
    #move i into the cluster or out
    new_edge = cluster_edges_dict[x]
    if x == cluster:
      global_edge_contr[0] = cluster_edges(i,x,get_cluster)
      new_edge += global_edge_contr[0]
      new_edge -= float(gamma*(math.comb(len(partition[x])+1,2)))
    else:
      new_edge -= cluster_edge_contribution[i][0]
      new_edge -= float(gamma*(math.comb(len(partition[x])-1,2)))
    #outfile.write(f"{x}: -{original} + {new_edge} ({cluster_edges_dict[x]}-{cluster_edge_contribution[i][0]})\n")
    change_sum += new_edge

  return change_sum

def cpm(gamma,partition):
  cpm_sum = 0
  for c in partition:
    cpm_sum += edges_in_cluster(partition,c) - gamma*(math.comb(len(partition[c]),2))
  return cpm_sum

#if u and v share an edge return 1
def a(u,v):
  if v in Neighbors[u]:
    return 1
  return 0

#given a node u, give the degree with respect to the cluster given
def cluster_edges(u, cluster, get_cluster):
  global global_edge_contr
  degree_sum = 0
  #Neightbors: given a key of a node, give its neighbors
  for neighboor in Neighbors[u]:
    #if a neighbor of u is in the cluster passed, increment
    if get_cluster[neighboor] == cluster:
      degree_sum += 1
      global_edge_contr[1].append(neighboor)
  return degree_sum

def move_nodes(h_old,g,partition,get_cluster,gamma,is_cpm):
  global cluster_degree_sum_dict
  global cluster_edge_contribution
  global cluster_edges_dict
  global exiting_modularity
  
  print("back in")
  start = process_time() 
  while True:
    total_change = 0
    for v in g.nodes():
      max_change = 0
      best_cluster = get_cluster[v] #put the original cluster as the best
      best_cluster_edges = cluster_edge_contribution[v]
      for cluster in partition:
        if len(partition[cluster]) > 0 and cluster != get_cluster[v] :
          change = change_in_cpm(v,cluster,gamma,partition,get_cluster) if (is_cpm) else change_in_modularity(v, cluster, gamma,get_cluster)
          best_cluster = cluster if (change > max_change) else best_cluster
          best_cluster_edges = global_edge_contr if (change > max_change) else best_cluster_edges
          max_change = change if (change > max_change) else max_change
      
      total_change += max_change
      #if max is bigger than 0 -> assign v to that cluster
      if max_change > 0:
        partition[get_cluster[v]].remove(v)  
        if (is_cpm):      
          cluster_edges_dict[get_cluster[v]] -= cluster_edge_contribution[v][0]
        else:
          cluster_degree_sum_dict[get_cluster[v]] -= len(Neighbors[v])
        
        for n in cluster_edge_contribution[v][1]:
          cluster_edge_contribution[n][0] -= 1
          cluster_edge_contribution[n][1].remove(v)
        
        get_cluster[v] = best_cluster
        partition[best_cluster].append(v)
        if (is_cpm):      
          cluster_edges_dict[get_cluster[v]] += best_cluster_edges[0]
        else:
          cluster_degree_sum_dict[best_cluster] += len(Neighbors[v])
        
        cluster_edge_contribution[v] = best_cluster_edges
        for n in best_cluster_edges[1]:
          cluster_edge_contribution[n][0] += 1
          cluster_edge_contribution[n][1].append(v)
    
    h_new = h_old + total_change
    if h_new <= h_old:
      print("gottem")
      #return new partition
      stop = process_time() 
      type_mod = "CPM" if (is_cpm) else "modularity"
      exiting_modularity = h_new
      output = f"louvain\t{start}\t{stop}\t{type_mod}\t{gamma}\t{h_new}\t{clean_partition(partition)}\t{disconnected_count(partition,get_cluster)}\t{badly_connected_count(g,partition,get_cluster,gamma,is_cpm)}" if (gamma == 1.0) else f"louvain\t{start}\t{stop}\t{type_mod}\t{gamma}\t{h_new}"
      outfile.write(output)
      outfile.write('\n')
      return partition
    else:
      h_old = h_new
      
def louvain(g,partition, get_cluster,gamma,is_cpm):
  threshold = 0.001
  h_old = cpm(gamma,partition) if (is_cpm) else modularity(gamma,partition)

  while True:
    partition = move_nodes(h_old,g,partition, get_cluster,gamma,is_cpm)
    if (exiting_modularity - h_old) < threshold:
      #return new partition
      return partition
    else:
      h_old = exiting_modularity

def move_nodes_fast(original,g,partition, get_cluster,gamma,is_cpm,print_a):
  global total_change_mnf
  global exiting_modularity
  global cluster_edge_contribution
  global cluster_degree_sum_dict
  #print(f"starting move_nodes_fast: partition:")
  start = process_time() 
  #loading the queue
  queue = [""] * len(g.nodes())
  in_queue = {}
  #visit nodes in a random order
  g_index = list(range(len(g.nodes())))
  random.shuffle(g_index)
  c = 0
  for i in g.nodes():
    queue[c] = i
    in_queue[i] = 1
    c += 1

  total_change = 0
  while len(queue) > 0:
    v = queue.pop()
    in_queue[v] = 0

    #determine the best community for node v
    max_change = 0
    best_cluster = get_cluster[v] #put the original cluster as the best
    best_cluster_edges = cluster_edge_contribution[v]
    for cluster in partition:
      if len(partition[cluster]) > 0 and cluster != get_cluster[v] :
        change = change_in_cpm(v,cluster,gamma,partition, get_cluster) if (is_cpm) else change_in_modularity(v, cluster, gamma, get_cluster)
        best_cluster = cluster if (change > max_change) else best_cluster
        best_cluster_edges = global_edge_contr if (change > max_change) else best_cluster_edges
        max_change = change if (change > max_change) else max_change

    total_change += max_change
    #if max is bigger than 0 -> assign v to that cluster
    if max_change > 0:
      for n in Neighbors[v]:
        #if a neighbor of v is not in C' and is not in the queue
        if get_cluster[n] != best_cluster and in_queue[n] == 0:
          queue.append(n)
        
      partition[get_cluster[v]].remove(v)  
      if (is_cpm):      
        cluster_edges_dict[get_cluster[v]] -= cluster_edge_contribution[v][0]
      else:
        cluster_degree_sum_dict[get_cluster[v]] -= len(Neighbors[v])
      
      for n in cluster_edge_contribution[v][1]:
        cluster_edge_contribution[n][0] -= 1
        cluster_edge_contribution[n][1].remove(v)
      
      get_cluster[v] = best_cluster
      partition[best_cluster].append(v)
      if (is_cpm):      
        cluster_edges_dict[get_cluster[v]] += best_cluster_edges[0]
      else:
        cluster_degree_sum_dict[best_cluster] += len(Neighbors[v])
      
      cluster_edge_contribution[v] = best_cluster_edges
      for n in best_cluster_edges[1]:
        cluster_edge_contribution[n][0] += 1
        cluster_edge_contribution[n][1].append(v)
  
  total_change_mnf = total_change
  #outputting needed data
  if print_a:
    stop = process_time() 
    h_new = original + total_change
    exiting_modularity = h_new
    type_mod = "CPM" if (is_cpm) else "modularity"
    output = f"leiden\t{start}\t{stop}\t{type_mod}\t{gamma}\t{h_new}\t{clean_partition(partition)}\t{0}\t{badly_connected_count(g,partition,get_cluster,gamma,is_cpm)}" if (gamma == 1.0) else f"leiden\t{start}\t{stop}\t{type_mod}\t{gamma}\t{h_new}"
    outfile.write(output)
    outfile.write('\n')
  #return new partition
  return partition #clean_partition(partition) 
   
#assign each node to it's cluster/community
def singleton_partition(g,get_cluster):
  partition = {}
  i = 0
  for node in g.nodes():
    #print(f"singleton partition chage {node} from {get_cluster[node]} to {i}")
    get_cluster[node] = i
    if i not in partition:
      partition[i] = []
    partition[i].append(node)
    i += 1
  return partition

def refine_partition(gamma,g,partition,get_cluster,is_cpm):
  global cluster_degree_sum_dict
  global cluster_edge_contribution
  global cluster_edges_dict
  p_refined = {}
  c_refined = {}
  index = 0
  count = 0
  for cluster in partition:
    cluster_singleton = {}
    index_2 = 0
    for node in partition[cluster]:
      if index_2 == 0:
        cluster_singleton[get_cluster[node]] = []
        cluster_singleton[get_cluster[node]].append(node)
      else:
        get_cluster[node] = max_partition_len + index
        cluster_singleton[get_cluster[node]] = []
        cluster_singleton[get_cluster[node]].append(node)
      index_2 += 1
      index += 1
      cluster_degree_sum_dict[get_cluster[node]] = len(Neighbors[node])
      cluster_edges_dict[get_cluster[node]] = 0
      cluster_edge_contribution[get_cluster[node]] = [0,[]]
    #come back here to handle the return here
    updated_data = merge_node_subset(gamma,g,cluster_singleton,get_cluster,cluster,is_cpm,partition)
    for community in updated_data:
      if len(updated_data[community]) > 0:
        p_refined[count] = updated_data[community]
        for v in updated_data[community]:
          c_refined[v] = count
      count += 1

  return [p_refined, c_refined]
    
def clean_partition(partition):
  count = 0
  for cluster in partition:
    if len(partition[cluster]) > 0:
      count += 1
  return count

def check_gamma_connected(gamma,community_p,cluster,original):
  s_minus_t_length = len(original[cluster]) - len(community_p)
  threshold =float(float(gamma) * (s_minus_t_length) * len(community_p))
  edge_num = 0
  for node in community_p:
    for neighbor in Neighbors[node]:
      if neighbor in original[cluster] and neighbor not in community_p:
        edge_num += 1
  return edge_num >= threshold

def merge_node_subset(gamma,g,partition,get_cluster,cluster,is_cpm,original):
  theta =0.1
  r = []
  #print(f"starting merge_node_subset of  with clust {cluster} ()")
  for v in partition:
    if check_gamma_connected(gamma,partition[v],cluster,original):
      r.append(v)
  for v in r:
    if len(partition[v]) == 1:
      elem = partition[v][0]
      T = []
      for community in partition:
        if check_gamma_connected(gamma,partition[community],cluster,original):
          T.append(community)
      probabilities = []
      for community in T:
        change = change_in_cpm(elem,community,gamma,partition,get_cluster) if (is_cpm) else change_in_modularity(elem, community, gamma,get_cluster)
        try:
          change = 0 if (change < 0) else math.exp((1/theta)*change) #exponential
        except:
          #at times I have gotten errors that the number within exp is too big 
          #so let take the max probability and increase by 10%, I know this isn't reliable, but I just want to get data
          change = max(probabilities)*1.10
        probabilities.append(change)
      #In Python, you can numpy.random.choice where these probabilities will be the final argument.
      random_cluster = int(np.random.choice(T,1,probabilities))
      partition[get_cluster[elem]].remove(elem)
      get_cluster[elem] = random_cluster #idk about this one
      partition[get_cluster[elem]].append(elem)
  return partition

def leiden(g,partition,get_cluster,gamma,is_cpm,a_print):
  threshold = 0.001
  if a_print:
    h_old = cpm(gamma,partition) if (is_cpm) else modularity(gamma,partition)
  else:
    h_old = exiting_modularity
  
  while True:
    print("going in move nodes fast")
    partition = move_nodes_fast(h_old,g,partition,get_cluster,gamma,is_cpm,a_print) #clean_partition(move_nodes_fast(g,partition,get_cluster,gamma,is_cpm))
    print("out of move nodes fast")
    #print(f"cleaned:")
    h_new = h_old + total_change_mnf
    if (h_new - h_old) < threshold:
      break
    else:
      h_old = h_new
      new_partition_cluster = refine_partition(gamma,g,partition, get_cluster,is_cpm) #find out how refine partition returns to partition
      partition = new_partition_cluster[0]
      get_cluster = new_partition_cluster[1]

def disconnected_count(partition, get_cluster):
  bad_sum = 0
  for cluster in partition:
    if len(partition[cluster]) > 0:
      node = partition[cluster][0]
      visited = {}
      visited[node] = 1
      queue = [node]
      for current_node in queue:
        queue.pop()
        for vertex in Neighbors[current_node]:
          if vertex not in visited and get_cluster[vertex] == cluster:
            visited[vertex] = 1
            queue.append(vertex)
      if len(visited) != len(partition[cluster]):
        bad_sum += 1
  return bad_sum

def badly_connected_count(g,partition,get_cluster,gamma,is_cpm):
  bad_sum = 0
  new_partition = partition
  new_get_cluster = get_cluster
  leiden(g,new_partition,new_get_cluster,gamma,is_cpm,False)
  for cluster in partition:
    if len(partition[cluster]) > 0:
      if not equal_partitions(new_partition[cluster],partition[cluster]):
        bad_sum += 1
  return bad_sum

def equal_partitions(l1, l2):
  #if lens are not equal, boom red flag
  if (len(l1) != len(l2)):
    return False

  l1.sort()
  l2.sort()
  for i in range(len(l1)):
    if (l1[i] != l2[i]):
        return False
  return True

def parse_human_network(k_count):
  global edges
  global Neighbors
  global Cluster_Dict
  global Cluster_List_Dict
  global max_partition_len
  global exiting_modularity
  global global_edge_contr
  global cluster_degree_sum_dict
  global cluster_edges_dict
  global cluster_edge_contribution
  #given a key of a node, its cluster is given
  Cluster_Dict = {}
  #given a key of a cluster, a list of the node within the cluster is given
  Cluster_List_Dict =  {}
  #given a key of a node, give its neighbors
  Neighbors = {}
  #num of edges in graph
  edges = 0
  #dict of cluster_degree_sum per cluster
  cluster_degree_sum_dict ={}
  #dict of edges within cluster
  cluster_edges_dict = {}
  #dict of edge contribution within cluster
  cluster_edge_contribution = {}
  #edges within a cluster
  global_edge_contr = [0, []]
  #modularity when we exit louvain/leiden
  exiting_modularity = 0
  #len of paritition
  max_partition_len = 0

  edge = {}
  g = nx.Graph()
  count = 0
  input_file = open('pathlinker-human-network.txt', 'r')
  lines = input_file.readlines()
  for line in lines[1:]:
    data = line.split('\t')
    if f"{data[0]}{data[1]}" not in edge and f"{data[1]}{data[0]}" not in edge:
      edge[f"{data[0]}{data[1]}"] = 1
      edge[f"{data[1]}{data[0]}"] = 1
      g.add_edge(data[0],data[1])
  
  g = nx.k_core(g, k=k_count, core_number=nx.core_number(g))
  count = 0
  for edge in g.edges():
    pair = 0
    first = ""
    second = ""
    for i in edge:
      edges += 1
      #record all nodes
      if i not in Cluster_Dict:
        Cluster_List_Dict[count] = [i]
        Cluster_Dict[i] = count
        cluster_edge_contribution[i] = [0,[]]
        Neighbors[i] = []
        count += 1
      if pair == 0:
        first = i
      else:
        second = i
      pair += 1
    Neighbors[first].append(second)
    Neighbors[second].append(first)
  max_partition_len = count
  input_file.close()
  print(len(g.nodes()))
  return g

def parse_string_network(k_count):
  global edges
  global Neighbors
  global Cluster_Dict
  global Cluster_List_Dict
  global max_partition_len
  global exiting_modularity
  global global_edge_contr
  global cluster_degree_sum_dict
  global cluster_edges_dict
  global cluster_edge_contribution
  #given a key of a node, its cluster is given
  Cluster_Dict = {}
  #given a key of a cluster, a list of the node within the cluster is given
  Cluster_List_Dict =  {}
  #given a key of a node, give its neighbors
  Neighbors = {}
  #num of edges in graph
  edges = 0
  #dict of cluster_degree_sum per cluster
  cluster_degree_sum_dict ={}
  #dict of edges within cluster
  cluster_edges_dict = {}
  #dict of edge contribution within cluster
  cluster_edge_contribution = {}
  #edges within a cluster
  global_edge_contr = [0, []]
  #modularity when we exit louvain/leiden
  exiting_modularity = 0
  #len of paritition
  max_partition_len = 0

  g = nx.Graph()
  input_file = open('STRING-network.csv', 'r')
  lines = input_file.readlines()
  for line in lines[1:]:
    data = line.split(',')
    #essentially record all edges -> getting neighbors
    g.add_edge(data[0],data[1])
    
  g = nx.k_core(g, k=k_count, core_number=nx.core_number(g))
  count = 0
  for edge in g.edges():
    pair = 0
    first = ""
    second = ""
    for i in edge:
      edges += 1
      #record all nodes
      if i not in Cluster_Dict:
        Cluster_List_Dict[count] = [i]
        Cluster_Dict[i] = count
        cluster_edge_contribution[i] = [0,[]]
        Neighbors[i] = []
        count += 1
      if pair == 0:
        first = i
      else:
        second = i
      pair += 1
    Neighbors[first].append(second)
    Neighbors[second].append(first)

  max_partition_len = count
  input_file.close()
  print(len(g.nodes()))
  return g

def data_collection():
  global Cluster_Dict
  global Cluster_List_Dict
  gammas = [0.1,1,10]
  for i in range(2):
    for q in range(2):
      for j in range(2):
        for gamma in gammas:
          g =  parse_human_network(15) if (j==0) else parse_string_network(15)
          name = "human" if (q == 0) else "string"
          outfile.write(f"{name}\n")
          outfile.write("algo\tstart\tstop\ttype\tgamma\tmod\tpartition_len\tdisconnected\tbadly_connected\n")
          if i == 0:
            print("louvain")
            louvain(g,Cluster_List_Dict,Cluster_Dict,gamma,(j==0))
          else:
            print("leiden")
            leiden(g,Cluster_List_Dict,Cluster_Dict,gamma,(j==0),True)
          print(f"finished {name} {i} {q} {j} {gamma}")

    outfile.write('\n')

def testing():
  global Cluster_Dict
  global Cluster_List_Dict
  global Neighbors
  global edges
   #given a key of a node, its cluster is given
  Cluster_Dict = {}
  #given a key of a cluster, a list of the node within the cluster is given
  Cluster_List_Dict =  {}
  #given a key of a node, give its neighbors
  Neighbors = {}
  #num of edges in graph
  edges = 0
  characts = ['a','b','c','d','e','f','g','h','i', 'j','k']
  Neighbors['a'] = ['d','e','b']
  Neighbors['b'] = ['a','c','f','h']
  Neighbors['c'] = ['b','k','f']
  Neighbors['d'] = ['a','g']
  Neighbors['e'] = ['a','g','h']
  Neighbors['f'] = ['b','c']
  Neighbors['g'] = ['d','e','h']
  Neighbors['h'] = ['g','e','i','j','b']
  Neighbors['i'] = ['h','j','k']
  Neighbors['j'] = ['h','i']
  Neighbors['k'] = ['c','i']

  for v in characts:
    cluster_edge_contribution[v] = [0,[]]
  
  g = nx.Graph()

  for i in range(len(characts)):
    edges += len(Neighbors[characts[i]])
    Cluster_Dict[characts[i]] = i
    if i not in Cluster_List_Dict:
      Cluster_List_Dict[i] = []
      g.add_node(characts[i])
    Cluster_List_Dict[i].append(characts[i])
  edges *= 0.5 #16, double checked this to be correct

  for v in g.nodes():
    for n in Neighbors[v]:
      g.add_edge(v,n)    
  
  louvain(g,Cluster_List_Dict,Cluster_Dict,0.1,True)
  print(cpm(0.1,Cluster_List_Dict))
  print(Cluster_List_Dict)

data_collection()
