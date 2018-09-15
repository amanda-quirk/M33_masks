using NearestNeighbors
using DelimitedFiles
using HDF5

rad2arcsec(r) = 3600 * rad2deg(r)

fileName = "/Users/amandaquirk/Documents/M33/Data/CFHT_data_for_isolation.hdf5"
dataFile = h5open(fileName, "r")
ra_all = read(dataFile, "RAs")
dec_all = read(dataFile, "Decs")
mag_all = read(dataFile, "mag")

# ra_all = ra_all[1:430]
# dec_all = dec_all[1:430]
# mag_all = mag_all[1:430]

# println( "Loading data")
# data =  readdlm( fileName )
nPoints = size(ra_all)[1] #indexing starts at 1!!
# # data = reshape( data, (2,nPoints))
# ra_all = data[:,1]
# dec_all = data[:,2]
# mag_all = data[:,3]

data_tree =  zeros( (2,nPoints)) #tree needs data to be in columns
for i in 1:nPoints
  data_tree[1,i] = ra_all[i]
  data_tree[2,i] = dec_all[i]
end

# Create tree
println( " Creating tree...")
kdtree = KDTree( data_tree; leafsize = 20 )
#
#
#
# ra = ra_all[1]
# dec = dec_all[2]
#
function get_angular_distance( point_1, point_2 )
  lambda_1, phi_1 = point_1*pi/180
  lambda_2, phi_2 = point_2*pi/180
  delta_lambda = lambda_1 - lambda_2
  dist = acos( sin(phi_1)*sin(phi_2) + cos( phi_1)*cos(phi_2)*cos(delta_lambda))
  return dist * 180 / pi * 3600
end
# phi_2, lambda_2

# point_1 = data[1,:]
# point_2 = data[2,:]
# dist_1 = get_angular_distance( point_1, point_2 )

const window_dist = 10 #make something a constant if you are going to use them in functions (faster)
rejected = zeros(nPoints)
function get_max_dist( i, nNeihg) #given the index and the number of neighbors to look for
  point = data_tree[:,i] #find neighbors for one point
  star_mag = mag_all[i]
  # println( point)
  # point_coord = ICRSCoords(point[1], point[2])  # inputs are ra, dec in radians
  idxs, dists = knn(kdtree, point, nNeihg, true) #true returns sorted by distance (not our distance)
  near_neigh = data_tree[:,idxs[2]] #first neighbor is itself so look for second neighbor
  near_dist = get_angular_distance( point, near_neigh)
  if near_dist > window_dist
    return
  end
  far_neigh = data_tree[:,idxs[end]] #end is last one in the array (farthest neighbor)
  far_dist = get_angular_distance(point, far_neigh)
  while far_dist < window_dist #probes all of the neighbors in our circle
    # println( " Searching more Neigbours")
    nNeihg *= 2 #increase the number of neighbors searched for to make sure we get all within the circle
    idxs, dists = knn(kdtree, point, nNeihg, true)
    far_neigh = data_tree[:,idxs[end]]
    far_dist = get_angular_distance(point, far_neigh)
  end
  #sum = 0
  for idx in idxs[2:end]
    point_neigh = data_tree[:,idx]
    dist = get_angular_distance( point, point_neigh)
    if dist > window_dist
      break
    end
    criteria = mag_all[idx] + (dist / 0.8)^2 - 3
    if criteria < star_mag #star NOT isolated if it has a neighbor that fufills this
        rejected[i] = 1
        return
    end

    # println( dist)
  end
  # max_dist_indx = idxs[end]
  # point_max_dist = data[max_dist_indx,:]

  # println( sum, sum )
  # dist_max = get_angular_distance( point, point_max_dist)
  # if dist < 10
    # println( dist )
  # end
  # max_dist_coord = ICRSCoords(ra_all[max_dist_indx], dec_all[max_dist_indx])
  # separ = separation.(point_coord, max_dist_coord)
  # separ = rad2deg()
  # return dist
end

println( " Computing distances...")

function find_distances( )
  nNeihg_base = 50
  for i in 1:nPoints
    if i % 5000 == 0
        println(i, " ", nPoints)
    end
    max_dist = get_max_dist(i, nNeihg_base)
  end
end
#
find_distances()
# println( time )

#writedlm( "/Users/amandaquirk/Documents/M33/Data/isolation_tag_GAIA_region1.txt", rejected)

fileName = "/Users/amandaquirk/Documents/M33/Data/CFHT_strict_isolated.hdf5"
dataFile = h5open(fileName, "w")
dataFile["isolation_tag"] = rejected
close(dataFile)
