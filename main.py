#! /usr/bin/python3
import svgpathtools
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches
from scipy.optimize import minimize
import numpy as np
import pdb

def complex2tuple(val):
	return (val.real, val.imag)

def convert_closed_path(svg_path):
	verts = [complex2tuple(svg_path[0].start)]
	print(svg_path[0].start.imag)
	codes = [Path.MOVETO]
	for segment in svg_path:
		if type(segment) is svgpathtools.CubicBezier:
			verts += [
				complex2tuple(segment.control1),
				complex2tuple(segment.control2),
				complex2tuple(segment.end),
			]
			codes += [Path.CURVE4]*3
			codes += [Path.MOVETO, Path.LINETO, Path.MOVETO,Path.LINETO ]
			verts += [complex2tuple(segment.start), complex2tuple(segment.control1), complex2tuple(segment.control2), complex2tuple(segment.end)]
			
		elif type(segment) is svgpathtools.Line:
			verts += complex2tuple(segment.end),
			codes += [Path.LINETO]
		else:
			print("Sorry. unexpected segment type")
			print(segment)
			exit(-1)
	return Path(verts, codes)

def tangent_paths(svg_path):
	verts = []
	codes = []
	seg_len = 20 #
	for segment in svg_path:
		if type(segment) is svgpathtools.CubicBezier:
			uv = segment.unit_tangent(1)
			tangent_end = segment.end + seg_len * uv
			uv = segment.unit_tangent(0)
			tangent_start = segment.start - seg_len * uv
			segment.point(0.5)
			verts += [
				complex2tuple(segment.end),
				complex2tuple(tangent_end),
				complex2tuple(segment.start),
				complex2tuple(tangent_start),	
			]
			codes += [Path.MOVETO, Path.LINETO, Path.MOVETO, Path.LINETO]
	return Path(verts, codes)

def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)

def angle_between(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'::

            >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

class Turn():
	def __init__(self, turn_type, start, end, arc, verts, profil_steps):
		self.turn_type = turn_type
		self.start = start
		self.end = end
		self.arc = arc
		self.verts = verts
		self.profil_steps = profil_steps

def calc_turn(start, vec, arc, r1, r2, turn_type, **kwargs):
	
	verts = []
	if r1 < 0.0 or r2 < 0.0:
		print("Warning: negative value(s) for radius r1 or r2") #TODO: raise exception
		return [start]
	# steps:
	total_length = (r1+r2)/2.0*arc # only approx. if r1 != r2
	N = 0
	if 'profil_step_length' in kwargs:
		N = int(total_length / kwargs.get('profil_step_length')) +1# number of steps
	else:
		N = kwargs.get('profil_steps')
	#print(N)
	#print(N)
	dRad = r2-r1	
	if N != 1:
		dRad = (r2-r1)/(N-1)
	dArcLength = total_length / N;
	if r1 != r2 and N != 1:
		# compensate for radius change
		dRad = (r2-r1)/(N-1) # TODO: N-1 seems a bit strange, keep an eye on it
		tmpAngle = 0.0
		tmpRadius = r1
		for n in range(N):
			tmpAngle += dArcLength / tmpRadius
			tmpRadius += dRad
		dArcLength = dArcLength * arc / tmpAngle # correct average length per step
	heading = np.arctan2(vec[1], vec[0])
	
	
	sign = 1
	if turn_type == 'lft':
		sign = -1
	ang_normal = heading - sign*np.pi/2.0
	newx, newy = start

	for n in range(0,N):
		seg_radius = r1+n*dRad
		seg_arc = dArcLength / seg_radius
		
		ceny = newy + seg_radius * np.sin(ang_normal)
		cenx =  newx + seg_radius * np.cos(ang_normal)
		#print([seg_radius, seg_arc, cenx, ceny, ang_normal])
		heading -= sign*seg_arc
		ang_normal -= sign*seg_arc
		newx = cenx - seg_radius * np.cos(ang_normal)
		newy = ceny - seg_radius * np.sin(ang_normal)
		if kwargs.get('all_verts', False):
			verts.append([newx, newy])
		#print(seg_radius)
	
	return Turn(turn_type, start, [newx, newy], arc, verts, N)

def err_dist(r, start, end_target, vec, arc, profil_step_length, turn_type):
	end_p = calc_turn_end(start, vec, arc, r[0], r[1], profil_step_length, turn_type, False)
	return np.linalg.norm(np.array(end_target)-np.array(end_p[0]))

def err_dist_jac(r, start, end_target, vec, arc, profil_step_length, turn_type):
	end_p0 = calc_turn_end(start, vec, arc, r[0]+2, r[1], profil_step_length, turn_type, False)
	end_p1 = calc_turn_end(start, vec, arc, r[0], r[1]+2, profil_step_length, turn_type, False)
	err_dist0 =  np.linalg.norm(np.array(end_target)-np.array(end_p0[0]))
	err_dist1 =  np.linalg.norm(np.array(end_target)-np.array(end_p1[0]))
	return np.array([err_dist0, err_dist1])

def neg_segment_curvature(x,segment):
	x = np.clip(x, 0.0, 1.0)
	print(x)
	return -segment.curvature(x[0])
def neg_segment_jacobian(x, segment):
	return -segment.derivative(x[0], n=1)


def cubic_bezier_inflections(cb):
	# taken from: https://pomax.github.io/bezierinfo/#inflections
	x0,y0 = complex2tuple(cb.start)
	x1,y1 = complex2tuple(cb.control1)
	x2,y2 = complex2tuple(cb.control2)
	x3,y3 = complex2tuple(cb.end)
	u = x3-3*x2+3*x1-x0
	v = x2-2*x1+x0
	w = x1-x0
	m = y3-3*y2+3*y1-y0
	n = y2-2*y1+y0
	o = y1-y0
	z = (n*w-o*v)
	y = (m*w-o*u)
	x = (m*v-n*u)
	if y**2-4*x*z < 0.0:
		return []
	t1 = -(np.sqrt(y**2-4*x*z)+y)/x/2.0E+0
	t2 = (np.sqrt(y**2-4*x*z)-y)/x/2.0E+0
	print([t1,t2])
	#t1 = -y + np.sqrt((y**2 - 4*x*z) / (2*x))
	#t2 = -y - np.sqrt((y**2 - 4*x*z) / (2*x))
	res = []
	if t1 >= 0.0 and t1 <= 1.0:
		res.append(t1)
	if t2 >= 0.0 and t2 <= 1.0:
		res.append(t2)
	print("Result:"+str(res))
	# align cubic bezier
	'''	rot_ang = np.arctan2(end[1], end[0])
	sin = np.sin(rot_ang)
	cos = np.cos(rot_ang)
	R = np.array([[cos, sin],[-sin, cos]])
	c1 = R.dot(c1)
	c2 = R.dot(c2)
	end = R.dot(end)
	pdb.set_trace()
	verts = [start, c1, c2, end]
	codes = [Path.MOVETO] + [Path.CURVE4]*3
	path = Path(verts, codes)
	patch = patches.PathPatch(path, facecolor='none', lw=1)
	fig = plt.figure()
	plt.gca().invert_yaxis()
	ax = fig.add_subplot(111)
	ax.add_patch(patch)

	ax.autoscale(enable=True, axis='both', tight=None)
	plt.axis('equal')
	plt.show()
	#pdb.set_trace()
	x2,y2 = (c1[0],c1[1])
	x3,y3 = (c2[0],c2[1])
	x4,y4 = (end[0],end[1])
	#print([cb.control1, cb.control2, cb.end])
	a = x3 * y2
	b = x4 * y2
	c = x2 * y3
	d = x4 * y3 
	x = 18*(-3*a + 2*b + 3*c - d)
	y = 18*(3*a - b - 3*c)
	z = 18*(c - a)
	#pdb.set_trace()
	print([x, y,z])
	if x == 0 and y != 0:
		print("Warn x == 0 and y != 0")
		return [-y]
	if (y**2 - 4*x*z) / (2*x) < 0.0:
		print("val < 0.0")
		return []
	t1 = -y + np.sqrt((y**2 - 4*x*z) / (2*x))
	t2 = -y - np.sqrt((y**2 - 4*x*z) / (2*x))
	print([t1,t2])
	res = []
	if t1 >= 0.0 and t1 <= 1.0:
		res.append(t1)
	if t2 >= 0.0 and t2 <= 1.0:
		res.append(t2)
	print("Result:"+str(res))'''
	return res
def draw_turn_to_figure(turn, fig, **kwargs):
	verts = [turn.start] + turn.verts
	codes = [Path.MOVETO]+[Path.LINETO]*(len(verts)-1)
	patch = patches.PathPatch(Path(verts,codes), facecolor='none', lw=1, edgecolor=kwargs.get('color', 'blue'))
	ax = fig.gca()
	ax.add_patch(patch)
	ax.autoscale(enable=True, axis='both', tight=None)
	plt.axis('equal')
	fig.show()

class MinimizerResult():
	def __init__(self, success, n_iter, err, res):
		self.success = success
		self.n_iter = n_iter
		self.error = err
		self.result = res

def custom_minimize(r, fig, **kwargs):
	#args=(start,end, vec, arc, profil_step_length, turn_type), options={'xtol': 1e-8, 'disp': True}):
	start, end_target, vec, arc, profil_step_length, turn_type = kwargs["args"]
	end_target = np.array(end_target)
	
	err_dist = 1e12 # dummy value
	r_new = r
	r_prev = r
	min_iter = 0
	max_iter = 1000
	n_iter = 0
	#pdb.set_trace()
	adjust_try_factor = np.array([1.0, 1.0])
	turn = None
	turn_opts = {'profil_step_length': profil_step_length}
	ftol = 1e-10
	while (n_iter < min_iter or err_dist > ftol) and n_iter < max_iter:

		turn = calc_turn(start, vec, arc, r_new[0], r_new[1], turn_type,
			all_verts=True, **turn_opts)
 		
		end_p = np.array(turn.end)
		err_dist_prev = err_dist
		err_dist = np.linalg.norm(end_target-end_p)
		
		#draw_turn_to_figure(turn, fig)
		####################
		adjust_try = err_dist *adjust_try_factor + np.array(r_new) * min(err_dist,0.1)
		# 
		turn0 = calc_turn(start, vec, arc, r_new[0]+adjust_try[0], r_new[1],
			turn_type, all_verts = False, **turn_opts) # other option: profil_steps = turn.profil_steps) #
		
		turn1 = calc_turn(start, vec, arc, r_new[0], r_new[1]+adjust_try[1],
			turn_type, all_verts = False, **turn_opts) # other option: profil_steps = turn.profil_steps) #

		end_p0 = np.array(turn0.end)
		end_p1 = np.array(turn1.end)
		v0 = end_p - end_p0
		dist0 = np.linalg.norm(v0)
		v1 = end_p - end_p1
		dist1 = np.linalg.norm(v1)
		uv0 = v0 / dist0
		uv1 = v1 / dist1
		a = np.eye(2, dtype=np.float64)
		a[:,0] = uv0
		a[:,1] = uv1
		r_prev = r_new
		b = end_p - end_target
		try:
			x = np.linalg.solve(a, b)
			if not np.allclose(np.dot(a, x), b):
				print("No solution found")
				adjust_try_factor[0] *=-1
				adjust_try_factor[1] *=-1
			else:
				r_new = 0.9*(x * adjust_try * [1/dist0, 1/dist1]) + r_new
		except  np.linalg.LinAlgError as err:
			adjust_try_factor[1] *=2
			adjust_try_factor[0] *=2
			#pdb.set_trace()
		if r_new[0] < 0.0 or r_new[1] < 0.0:
			# if we attempt to use a neg. radius, try opposite turn type 
			#if turn_type == 'rgt':
			#	turn_type = 'lft'
			#else:
			#	turn_type = 'rgt'
			#print("Switching turn type to '%s'"%(turn_type))
			print("Attempting to use neg radius :(")
			if r_new[0] < 0.0:
				r_new[0] = profil_step_length #1.0 #0.9*r_prev[0]
				#adjust_try_factor[0] *=-1
				#r_new[1] = 0.9*r_new[1]
			if r_new[1] < 0.0:
				#adjust_try_factor[1] *=-1
				#r_new = r_prev
				r_new[1] = profil_step_length # 0.9*r_prev[1]
		if err_dist_prev < err_dist and 'profil_step_length' in turn_opts:
			# not converging
			# clamp the number of steps, otherwise minimizer might oscillate between 2 solutions
			# if the number of steps changes between 2 optimisation steps
			min_iter = n_iter + 10
			turn_opts.pop('profil_step_length')
			turn_opts['profil_steps'] = turn.profil_steps
		if n_iter >= min_iter and err_dist > profil_step_length and 'profil_steps' in turn_opts:
			# remove the clamping for the number of steps again
			turn_opts.pop('profil_steps')
			turn_opts['profil_step_length'] = profil_step_length
		'''if np.allclose(r_new, r_prev, atol=profil_step_length/4) and 'profil_step_length' in turn_opts:
			# clamp the number of steps, otherwise minimizer might oscillate between 2 solutions
			# if the number of steps changes between 2 optimisation steps
			print("Clamp profil steps")
			turn_opts.pop('profil_step_length')
			turn_opts['profil_steps'] = turn.profil_steps
				#r_new[0] = 0.9*r_new[0]
			# restart with positive radius
			#chord = np.linalg.norm(np.array(end_target)-np.array(start))
			#r_new[0] = r_new[1] = chord / (2.0* np.sin(arc/2.0))'''
		#if  np.allclose(end_target[0], 361.19357715):
		#	pdb.set_trace()
		n_iter +=1
		print("Err dist: "+str(err_dist)+" Radius: "+str(r_new))
		#pdb.set_trace()
	if 'profil_steps' in turn_opts:
		print("n_iter ="+str(n_iter)+" max_iter ="+str(max_iter))
		print("number of steps clamped to "+str(turn_opts['profil_steps']))
	success = err_dist < ftol
	
		
	return MinimizerResult(success, n_iter, err_dist, turn)

def angle_unwrap_2pi(ang):
	while ang < 0:
		ang += 2*np.pi
	return ang
		
def cubic_bezier_segment_to_turn(segment, t1, t2, verts, codes, fig):
	uv1 = complex2tuple(segment.unit_tangent(t2))
	uv0 = complex2tuple(segment.unit_tangent(t1))
	arc = angle_between(uv0, uv1)
	
	uv1_ang = angle_unwrap_2pi(np.arctan2(uv1[1], uv1[0]))
	uv0_ang = angle_unwrap_2pi(np.arctan2(uv0[1], uv0[0]))
	print("Start/end tangents: "+str([uv0_ang/np.pi*180.0, uv1_ang/np.pi*180.0 ]))
	# cw or ccw:
	turn_dir = uv0[0] * uv1[1] - uv0[1] * uv1[0]
	turn_type = 'rgt'
	if turn_dir > 0.0:
		turn_type = 'lft'
	print("-> Turn type: %s"%(turn_type))
	start = complex2tuple(segment.point(t1))
	end = complex2tuple(segment.point(t2))
	vec = uv0
	# initial guess:
	# assume r1 and r2 are equal, use the formula for a 'chord of a circle' 
	chord = np.linalg.norm(np.array(end)-np.array(start))
	r1 = r2 = chord / (2.0* np.sin(arc/2.0))
	profil_step_length = 4.0
	#res = minimize(err_dist, [r1,r2], args=(start,end, vec, arc, profil_step_length, turn_type), method='Newton-CG', options={'xtol': 1e-8, 'disp': False, 'maxiter':10000, 'fatol':1e-10, 'xatol':1e-10}, tol=1e-10, jac=err_dist_jac)
	#turn_verts = calc_turn_end(start, vec, arc, res.x[0], res.x[1], profil_step_length, turn_type, True)
	res = custom_minimize([r1,r2], fig, args=(start,end, vec, arc, profil_step_length, turn_type))
	#print(r_new)
	if not res.success:
		print("NO SOLUTION")
		draw_turn_to_figure(res.result, fig, color='red')
	else:
		draw_turn_to_figure(res.result, fig, color='green')
	return res	
	

def fit_torcs_segments(svg_path, fig):
	verts = []
	codes = []
	n_iter = []
	for n,segment in enumerate(svg_path):
		if n < 1 or n > 2:
			#continue
			pass
		if type(segment) is svgpathtools.CubicBezier:
			#max_curvature = minimize(neg_segment_curvature, [0.5], args=(segment,), method='L-BFGS-B', options={'xtol': 1e-5, 'disp': True}, bounds=((0,1.0),))
			#print(max_curvature.x[0])
			#t_end = 0.6#max_curvature.x[0]
			profil_step_length = 4.0
			t_vals = [0.0] + cubic_bezier_inflections(segment) +[1.0]
			i = 0
			n_insert_inbetween = 0 # number of new
			while i < len(t_vals)-1:
				t_start = t_vals[i]
				t_end = t_vals[i+1]
				
				t_dist = (t_end - t_start)/(n_insert_inbetween+1.0)
				i += 1
				for j in range(1,n_insert_inbetween +1):
					t_vals.insert(i, t_start + t_dist*j)
					i +=1
					#pdb.set_trace()
				
			print(t_vals)
			i = 0
			while i < len(t_vals)-2:
				t_start = t_vals[i]
				t_end = t_vals[i+1]
				if segment.length(t_start  , t_end) < profil_step_length:
					print("Warning: removing turn segment (length < profil_step_length)")
					t_vals.pop(i+1)
				else:
					i += 1
					
			print("t values: "+str(t_vals))
			for i in range(len(t_vals)-1):
				t_start = t_vals[i]
				t_end = t_vals[i+1]
				ax = fig.gca()
				end_p = complex2tuple(segment.point(t_end))
				ax.plot(end_p[0], end_p[1], 'x')
				res = cubic_bezier_segment_to_turn(segment, t_start, t_end, verts, codes, fig)
				n_iter.append(res.n_iter)
	print(np.mean(n_iter))
	#print(len(verts))
	#print(len(codes))
	return #Path(verts, codes)


paths, attributes = svgpathtools.svg2paths('track_layout.svg')


fig = plt.figure()
plt.gca().invert_yaxis()
ax = fig.add_subplot(111)

for svg_path in paths[1:2]:
	mat_path = convert_closed_path(svg_path)
	patch = patches.PathPatch(mat_path, facecolor='none', lw=1)
	ax.add_patch(patch)

for svg_path in paths[1:2]:
	mat_path = tangent_paths(svg_path)
	patch = patches.PathPatch(mat_path, facecolor='none', lw=1, edgecolor='red')
	ax.add_patch(patch)

for svg_path in paths[1:2]:
	torcs_path = fit_torcs_segments(svg_path, fig)
	#patch = patches.PathPatch(torcs_path, facecolor='none', lw=1, edgecolor='green')
	#ax.add_patch(patch)
ax.autoscale(enable=True, axis='both', tight=None)
plt.axis('equal')
plt.show()
