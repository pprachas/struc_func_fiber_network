import meshio
import numpy as np
import matplotlib.pyplot as plt
#----------Import functions---------------------#
import sys
sys.path.append('../../utils')
import meshing 
import graph


#-------------Parameters-------------------#
W = 10000 # window width
H = 10000 # window height
L = W*10 # length 
n = 60 #number of fibers
#-------------Set seed---------------------#
rng = np.random.RandomState(0)

#------Get Line centers and orientation----#
centers_x = rng.uniform(0,W,n)
centers_y = rng.uniform(0,H,n)
theta = rng.uniform(0,2*np.pi,n)
#--------------Create initial lines--------#
end_1 = np.array([centers_x + (L/2)*np.cos(theta),centers_y + (L/2)*np.sin(theta)])
end_2 =np.array([centers_x + (L/2)*np.cos(theta-np.pi),centers_y + (L/2)*np.sin(theta-np.pi)])
#-----------Trim fibers outside window-----#

# x direction
left_bound_1 = np.where(end_1[0,:]<0)
for ii in left_bound_1:
    end_1[1,ii] = centers_y[ii] + (-centers_x[ii])*((end_1[1,ii]-centers_y[ii])/(end_1[0,ii]-centers_x[ii]))
    end_1[0,ii] = 0.0

right_bound_1 = np.where(end_1[0,:]>W)
for ii in right_bound_1:
    end_1[1,ii] = centers_y[ii] + (W-centers_x[ii])*((end_1[1,ii]-centers_y[ii])/(end_1[0,ii]-centers_x[ii]))
    end_1[0,ii] = W

left_bound_2 = np.where(end_2[0,:]<0)
for ii in left_bound_2:
    end_2[1,ii] = centers_y[ii] + (-centers_x[ii])*((end_2[1,ii]-centers_y[ii])/(end_2[0,ii]-centers_x[ii]))
    end_2[0,ii] = 0.0
    
right_bound_2 = np.where(end_2[0,:]>W)
for ii in right_bound_2:
    end_2[1,ii] = centers_y[ii] + (W-centers_x[ii])*((end_2[1,ii]-centers_y[ii])/(end_2[0,ii]-centers_x[ii]))
    end_2[0,ii] = W

# y direction

lower_bound_1 = np.where(end_1[1,:]<0)
for ii in lower_bound_1:
    end_1[0,ii] = centers_x[ii] + (-centers_y[ii])*((end_1[0,ii]-centers_x[ii])/(end_1[1,ii]-centers_y[ii]))
    end_1[1,ii] = 0.0

upper_bound_1 = np.where(end_1[1,:]>H)
for ii in upper_bound_1:
    end_1[0,ii] = centers_x[ii] + (H-centers_y[ii])*((end_1[0,ii]-centers_x[ii])/(end_1[1,ii]-centers_y[ii]))
    end_1[1,ii] = H

lower_bound_2 = np.where(end_2[1,:]<0)
for ii in lower_bound_2:
    end_2[0,ii] = centers_x[ii] + (-centers_y[ii])*((end_2[0,ii]-centers_x[ii])/(end_2[1,ii]-centers_y[ii]))
    end_2[1,ii] = 0.0

upper_bound_2 = np.where(end_2[1,:]>H)  
for ii in upper_bound_2:
    end_2[0,ii] = centers_x[ii] + (H-centers_y[ii])*((end_2[0,ii]-centers_x[ii])/(end_2[1,ii]-centers_y[ii]))
    end_2[1,ii] = H

#----------------Find intersection points--------------------------#
intersect = []
line_num = []
for ii in range(0,n):
    line1 = np.array([end_1[1,ii]-end_2[1,ii],end_2[0,ii]-end_1[0,ii], end_2[0,ii]*end_1[1,ii] - end_1[0,ii]*end_2[1,ii]]) # coefficients

    for jj in range(0,n):
        if ii != jj: # cramer's rule
            line2 = np.array([end_1[1,jj]-end_2[1,jj],end_2[0,jj]-end_1[0,jj], end_2[0,jj]*end_1[1,jj] - end_1[0,jj]*end_2[1,jj]])
            
            D  = line1[0] * line2[1] - line1[1] * line2[0]
            Dx = line1[2] * line2[1] - line1[1] * line2[2]
            Dy = line1[0] * line2[2] - line1[2] * line2[0]
            if D != 0:
                x = Dx / D
                y = Dy / D
                
                # Apppend only points located on segment
                if np.minimum(end_1[0,jj],end_2[0,jj])<=x<=np.maximum(end_1[0,jj],end_2[0,jj]) and np.minimum(end_1[1,jj],end_2[1,jj])<=y<=np.maximum(end_1[1,jj],end_2[1,jj]):
                    if np.minimum(end_1[0,ii],end_2[0,ii])<=x<=np.maximum(end_1[0,ii],end_2[0,ii]) and np.minimum(end_1[1,ii],end_2[1,ii])<=y<=np.maximum(end_1[1,ii],end_2[1,ii]):
                        intersect.append([x,y])
                        line_num.append([ii,jj])

intersect = np.array(intersect)
#intersect, idx = np.unique(intersect, axis = 0, return_index = True)

line_num = np.array(line_num)

mesh_points_beams = []
mesh_points_mold = []
for ii in range(0,n):
    line_points = line_num[:,0]==ii
    points_on_line = np.vstack((end_1[:,ii], end_2[:,ii], intersect[line_points,:]))
    sort_idx = np.argsort(points_on_line[:,0])
    mesh_points_beams.append(points_on_line[sort_idx])
    # mesh_points_mold.append([end_1[:,ii], end_2[:,ii]])

dist = []
for ii in range(0,len(mesh_points_beams)):
    for jj in range(0,len(mesh_points_beams[ii])-1):
        plt.plot([mesh_points_beams[ii][jj,0],mesh_points_beams[ii][jj+1,0]],[mesh_points_beams[ii][jj,1],mesh_points_beams[ii][jj+1,1]])
        dist.append(np.sqrt((mesh_points_beams[ii][jj,0]-mesh_points_beams[ii][jj+1,0])**2+(mesh_points_beams[ii][jj,1]-mesh_points_beams[ii][jj+1,1])**2))

print(np.mean(np.array(dist)))
plt.axis('equal')
plt.show()
meshing.create_mesh('.','mikado',mesh_points_beams,char_length = L/50)

print(np.mean(np.array(dist)))
