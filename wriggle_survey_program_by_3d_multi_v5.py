#wriggle_survey_program_by_3d_multi_v5
import numpy as np
from scipy import optimize
import csv

# Convert Angle : แปลงมุม (Credit Prajuab Riabroy's Blog)
PI = np.pi
DEG2RAD = PI / 180.0
RAD2DEG = 180.0 / PI

################################## Function List ##################################

#Write csv file
def WriteCSV (export_best_fit_result):
	with open ('export_best_fit_result.csv', 'a', newline='', encoding='utf-8') as file:
		#fw is file writer
		fw = csv.writer(file)
		fw.writerow(export_best_fit_result)

#Direction Function (Distance and Azimuth)
def Direction(Estart, Nstart, Eend, Nend):
	dE = Eend - Estart
	dN = Nend - Nstart

	Dist = np.sqrt(dE**2 + dN**2)
	A = np.arctan(dE / dN) * RAD2DEG

	if dN < 0:
		Azi = 180 + A
	elif dE < 0:
		Azi = 360 + A
	else:
		Azi = A		
	return Dist, Azi		

#Pitching Function
def Pitch(Chstart, Elstart, Chend, Elend):
	P = (Elend - Elstart) / (Chend - Chstart)
	return P

#Horizontal Deviation Function
def DHz (Ed, Nd, Azd, Ea, Na):
	l = np.sqrt((Ea - Ed)**2 + (Na - Nd)**2)
	Aza = Direction(Ed, Nd, Ea, Na)
	AErr = (Aza[1] - Azd) * PI / 180.0	

	CH = l * np.cos(AErr)
	OS = l * np.sin(AErr)
	return CH, OS

#Vertical Deviation	Function
def DVt(Chd, Eld, P, Chf, Elf):
	El = Eld + P * (Chf - Chd)
	Vt = Elf - El  
	return Vt		

################################ Best-Fit Circle Center ###############################

#As built Coordinated Data by 3D
num = 8 #(1 ring / 8 points)
Data = np.loadtxt('1_import_WRS_R0toR917.csv', delimiter=',', skiprows=1, dtype=None) #Import data.csv

#Count the data in list (Total data)
for count in enumerate(Data):
	count
total_ring = int((count[0]+1)/num)

Data_spl = np.vsplit(Data, total_ring) #Split data

print('Best Fit Center Result')
WRS = [] #Best Fit Center append
for i in range(total_ring): 
	Data_spl[i]
	P = np.asarray(Data_spl[i][:,0]) #Column 1 is point (1 ring / 8 points)
	E = np.asarray(Data_spl[i][:,1]) #column 2 is Easting
	N = np.asarray(Data_spl[i][:,2]) #Column 3 is Northing
	Z = np.asarray(Data_spl[i][:,3]) #Column 4 is Elevation
	Ring_no = int(np.mean(P)) #Ring no. of data set

	#Linear Regression by Least square
	m, b = np.polyfit(E, N, 1)

	EMin = np.amin(E) * 0.999999
	EMax = np.amax(E) * 1.0000005

	#P and Q point on line fitting
	EP, NP = EMin, m * EMin + b
	EQ, NQ = EMax, m * EMax + b

	#Coordinates on PQ Line and Local coordinates xi, yi
	Local_xy = []
	for Pi, Ei, Ni, Zi in Data_spl[i]:
		E0 = (m * Ni + Ei - m * b) / (m**2 + 1)
		N0 = m * ((m * Ni + Ei - m * b) / (m**2 + 1)) + b

		x = np.sqrt((EP - E0)**2 + (NP - N0)**2)
		y = Zi + 100

		Local_xy.append([x, y])

	#Best-fit Circle by least square
	Local_xy = np.array(Local_xy) #Local Coordinates 2D (list)

	x = np.asarray(Local_xy[:,0])
	y = np.asarray(Local_xy[:,1])

	# coordinates of the barycenter
	x_m = np.mean(x)
	y_m = np.mean(y)

	def calc_R(xc, yc):
	    #calculate the distance of each 2D points from the center (xc, yc)
	    return np.sqrt((x-xc)**2 + (y-yc)**2)

	def f_2(c):
	    #calculate the algebraic distance between the data points and the mean circle centered at c=(xc, yc)
	    Ri = calc_R(*c)
	    return Ri - np.mean(Ri)

	center_estimate = x_m, y_m
	center_2, ier = optimize.leastsq(f_2, center_estimate)

	xc_2, yc_2 = center_2 #Xc, Yc are Center local coordinates
	Ri_2       = calc_R(*center_2) #Ri_2 is each radius
	Rc_2       = np.mean(Ri_2) #Rc_2 is averaged raduis
	sd = np.std(Ri_2, None) #Standard Deviation of Radius

	#Transform local coordinates center (Xc, Yc) to coordinates center (Ec, Nc, Zc)
	Q = np.arctan(m) #Q angle in radians

	Ec = EP + xc_2 * np.cos(Q)
	Nc = NP + xc_2 * np.sin(Q)
	Zc = yc_2 - 100

	WRS.append([Ring_no, Ec, Nc, Zc, Rc_2])

########################## Best-Fit Circle Center Deviation from Tunnel Axis ##########################

#Import Data Tunnel Axis(DTA)
DTA = np.loadtxt('2_import_DTA_SWLtoSST_dl.csv', delimiter=',', skiprows=1, dtype=None)

#Best fit center Data (List)
F = np.array(WRS) 

col_name = ['RING NO.', 'EASTING', 'NORTHING', 'ELEVATION', 'CHAINAGE', 'dHz', 'dVt', 'R_Average']
WriteCSV(col_name)
for RF, EF, NF, ZF, rF in F:
	#Index point (A, B, C) from Tunnel Axis
	l = [] #distance from i to F (Best Fit point)

	for Pt, Ct, Et, Nt, Z in DTA:	
		r = np.sqrt((Et - EF)**2 + (Nt - NF)**2)
		l.append([r])

	B_pt = l.index(min(l)) #minimum distance position in array l[]
	Index_A = DTA[B_pt - 1] #Point A[P, C, E, N, Z] 
	Index_B = DTA[B_pt] #Point B[P, C, E, N, Z] (nearly)
	Index_C = DTA[B_pt + 1] #Point C[P, C, E, N, Z]

	#Distance AF, CF Azimuth AB, BC and Pitching AB, BC
	DireAF = Direction(Index_A[2], Index_A[3], EF, NF)
	DireCF = Direction(Index_C[2], Index_C[3], EF, NF)
	DireAB = Direction(Index_A[2], Index_A[3], Index_B[2], Index_B[3])
	DireBC = Direction(Index_B[2], Index_B[3], Index_C[2], Index_C[3])
	PitAB = Pitch(Index_A[1], Index_A[4], Index_B[1], Index_B[4])
	PitBC = Pitch(Index_B[1], Index_B[4], Index_C[1], Index_C[4])

	#Deviation of F point (Best Fit) from Tunnel Axis
	if DireAF[0] < DireCF[0]:
		DevHz = DHz(Index_B[2], Index_B[3], DireAB[1], EF, NF)
		ChF = Index_B[1] + DevHz[0]
		HzF = DevHz[1] 
		VtF = DVt(Index_B[1], Index_B[4], PitAB, ChF, ZF)
	else:
		DevHz = DHz(Index_B[2], Index_B[3], DireBC[1], EF, NF)
		ChF = Index_B[1] + DevHz[0]
		HzF = DevHz[1] 
		VtF = DVt(Index_B[1], Index_B[4], PitBC, ChF, ZF)

	print('R{:.0f} E: {:.3f} N: {:.3f}, Z: {:.3f} Ch: {:.3f} Hz: {:.3f}, Vt: {:.3f}, R: {:.3f}'.format(RF, EF, NF, ZF, ChF, HzF, VtF, rF))

	Data_Result = [RF, EF, NF, ZF, ChF, HzF, VtF, rF]
	WriteCSV(Data_Result)