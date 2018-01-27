#! /usr/bin/python3
import svgpathtools
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches
from scipy.optimize import minimize
import numpy as np
import pdb
import torcs_trackxml
from track_segment import Turn
from timeit import default_timer as timer


   
def complex2tuple(val):
	return (-val.real, val.imag)

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



def calc_turn(start, heading, arc, r1, r2, turn_type, **kwargs):
	
	
	if r1 < 0.0 or r2 < 0.0:
		print("Warning: negative value(s) for radius r1 or r2") #TODO: raise exception
		return  Turn(turn_type, start, start, arc, verts, N)

	dtype = kwargs.get('dtype', np.float64)
	# convert all input to desired dtype
	start = np.array(start, dtype=dtype)
	heading = dtype(heading)
	arc = dtype(arc)
	r1 = dtype(r1)
	r2 = dtype(r2)

	

	# required number of steps
	total_length = (r1+r2)/dtype(2.0)*arc # only approx. if r1 != r2
	N = 0
	if 'profil_step_length' in kwargs:
		N = int(total_length / kwargs.get('profil_step_length')) +1# number of steps
	else:
		N = kwargs.get('profil_steps')

	# radius change per step
	dRad = r2-r1	
	if N != 1:
		dRad = dtype( (r2-r1)/(N-1) )

	dArcLength = dtype(total_length / N);
	
	if r1 != r2 and N != 1:
		# compensate for radius change
		#dRad = (r2-r1)/(N-1) # TODO: N-1 seems a bit strange, keep an eye on it
		tmpAngle = dtype(0.0)
		tmpRadius = r1
		for n in range(N):
			tmpAngle += dArcLength / tmpRadius
			tmpRadius += dRad
		dArcLength = dArcLength * arc / tmpAngle # correct average length per step

	
	
	sign = dtype(1)
	if turn_type == 'lft':
		sign = dtype(-1)


	ang_normal = heading - dtype(sign*np.pi/2.0)
	newx, newy = start
	verts = []
	for n in range(0,N):
		seg_radius = r1+dtype(n*dRad)
		seg_arc = dArcLength / seg_radius
		
		ceny = newy + seg_radius * np.sin(ang_normal)
		cenx =  newx + seg_radius * np.cos(ang_normal)
		ang_normal -= sign*seg_arc
		newx = cenx - seg_radius * np.cos(ang_normal)
		newy = ceny - seg_radius * np.sin(ang_normal)
		if kwargs.get('all_verts', False):
			verts.append([newx, newy])
	if dtype is np.float32:
		print(isinstance(dArcLength, np.float32))
		print(isinstance(ang_normal, np.float32))
		print(isinstance(total_length, np.float32))
		#pdb.set_trace()
		
	return Turn(turn_type, start, [newx, newy], arc, verts, N)

def err_dist(r, start, end_target, vec, arc, profil_step_length, turn_type):
	end_p = calc_turn_end(start, vec, arc, r[0], r[1], profil_step_length, turn_type, False)
	return np.linalg.norm(np.array(end_target)-np.array(end_p[0]))


def cubic_bezier_inflections(cb):
	# see the wx-maxima sheet in material for details
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
	res = []
	if t1 >= 0.0 and t1 <= 1.0:
		res.append(t1)
	if t2 >= 0.0 and t2 <= 1.0:
		res.append(t2)
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
	start, end_target, heading, arc, profil_step_length, turn_type = kwargs["args"]
	end_target = np.array(end_target, dtype=np.float32)
	
	err_dist = 1e12 # high start value
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

		turn = calc_turn(start, heading, arc, r_new[0], r_new[1], turn_type,
			all_verts=True, **turn_opts)
 		
		end_p = np.array(turn.end)
		err_dist_prev = err_dist
		err_dist = np.linalg.norm(end_target-end_p)
		# for debugging:
		#draw_turn_to_figure(turn, fig)
		####################
		adjust_try = err_dist *adjust_try_factor + np.array(r_new) * min(err_dist,0.1)
		# 
		turn0 = calc_turn(start, heading, arc, r_new[0]+adjust_try[0], r_new[1],
			turn_type, all_verts = False, **turn_opts) # other option: profil_steps = turn.profil_steps) #
		
		turn1 = calc_turn(start, heading, arc, r_new[0], r_new[1]+adjust_try[1],
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
			print("Attempting to use neg radius :(")
			if r_new[0] < 0.0:
				r_new[0] = profil_step_length
			if r_new[1] < 0.0:
				r_new[1] = profil_step_length
		if err_dist_prev < err_dist and 'profil_step_length' in turn_opts:
			# not converging
			# clamp the number of steps, otherwise minimizer might oscillate between 2 solutions 
			# ... if the number of steps changes between 2 optimisation steps
			min_iter = n_iter + 10
			turn_opts.pop('profil_step_length')
			turn_opts['profil_steps'] = turn.profil_steps
		if n_iter >= min_iter and err_dist > profil_step_length and 'profil_steps' in turn_opts:
			# remove the clamping for the number of steps again
			turn_opts.pop('profil_steps')
			turn_opts['profil_step_length'] = profil_step_length
		n_iter +=1
		print("Err dist: "+str(err_dist)+" Radius: "+str(r_new))
		#pdb.set_trace()
	if 'profil_steps' in turn_opts:
		print("n_iter ="+str(n_iter)+" max_iter ="+str(max_iter))
		print("number of steps clamped to "+str(turn_opts['profil_steps']))

	success = err_dist < ftol
	turn.set_radius(*r_new)	
	return MinimizerResult(success, n_iter, err_dist, turn)
		
def cubic_bezier_segment_to_turn(segment, t1, t2, fig):
	uv0 = np.array(complex2tuple(segment.unit_tangent(t1)), dtype=np.float32)
	uv1 = np.array(complex2tuple(segment.unit_tangent(t2)), dtype=np.float32)
	arc = angle_between(uv0, uv1)
	
	# figure out whether turn is left or right, cw or ccw:
	turn_dir = uv0[0] * uv1[1] - uv0[1] * uv1[0]
	turn_type = 'rgt'
	if turn_dir > 0.0:
		turn_type = 'lft'

	print("-> Turn type: %s"%(turn_type))
	start = complex2tuple(segment.point(t1))
	end = complex2tuple(segment.point(t2))

	heading = np.arctan2(uv0[1], uv0[0])
	# initial guess: assume r1 and r2 are equal, use the formula for a 'chord of a circle' 
	chord = np.linalg.norm(np.array(end)-np.array(start))
	r1 = r2 = chord / (2.0* np.sin(arc/2.0))

	# default profile step length:
	profil_step_length = 4.0

	# alternative option: use the scipy minimizer
	#res = minimize(err_dist, [r1,r2], args=(start,end, vec, arc, profil_step_length, turn_type), method='Newton-CG', options={'xtol': 1e-8, 'disp': False, 'maxiter':10000, 'fatol':1e-10, 'xatol':1e-10}, tol=1e-10, jac=err_dist_jac)
	#turn_verts = calc_turn_end(start, vec, arc, res.x[0], res.x[1], profil_step_length, turn_type, True)
	res = custom_minimize([r1,r2], fig, args=(start,end, heading, arc, profil_step_length, turn_type))

	return res	
	

def fit_torcs_segments(svg_path, fig):
	segments = []
	n_iter = []
	for segment in svg_path:
		if type(segment) is svgpathtools.CubicBezier:
			profil_step_length = 4.0
			t_vals = [0.0] + cubic_bezier_inflections(segment) +[1.0]
			i = 0
			n_insert_inbetween = 1 # number of new
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
				#ax = fig.gca()
				end_p = complex2tuple(segment.point(t_end))
				#ax.plot(end_p[0], end_p[1], 'x')
				res = cubic_bezier_segment_to_turn(segment, t_start, t_end, fig)
				color = 'green'
				if not res.success:
					color = 'red'
				#draw_turn_to_figure(res.result, fig, color=color )
				n_iter.append(res.n_iter)
				segments.append(res.result)
	print(np.mean(n_iter))
	return segments

def fix_closing_for_float32(segments):
	fixed_segments = []
	profil_step_length = 4.0
	last_turn = None
	last_turn_index = None
	last_turn_target_end =   np.array((0,0), dtype=np.float32)
	for i, segment in enumerate(reversed(segments)):
		if not segment.type == 'str':
			last_turn = segment
			last_turn_index = len(segments)-1-i
			break
		else:
			last_turn_target_end[0] -= segment.len
	# calculate last turn start position with float32 precision
	pos = np.array((0,0), dtype=np.float32)
	heading =np.float32(0.0)
	for i, segment in  enumerate(segments):
		if i == last_turn_index:
			break
		turn = calc_turn(pos, heading, segment.arc, segment.radius, segment.end_radius, segment.type, profil_steps=segment.profil_steps, dtype=np.float32)
		pos = turn.end
		print(turn.end)
		print( np.array(turn.end, dtype=np.float32))
		if turn.type == 'lft':
			heading += np.float32(turn.arc)
		else:
			heading -= np.float32(turn.arc)

	old_last_turn = calc_turn(pos, heading, last_turn.arc, last_turn.radius, last_turn.end_radius, last_turn.type, profil_steps=last_turn.profil_steps, dtype=np.float32)
	print("float32 drift: "+str(np.linalg.norm(np.array(old_last_turn.end)-np.array(last_turn.end))))
	print(pos)
	pdb.set_trace()
	res = custom_minimize([last_turn.radius, last_turn.end_radius], None, args=(pos, last_turn_target_end, heading, last_turn.arc, profil_step_length, last_turn.type))
	if res.success:
		segments.remove(last_turn)
		segments.insert(last_turn_index, res.result)
		print(res.result.end)
		turn32 = calc_turn(pos, heading, res.result.arc, res.result.radius, res.result.end_radius, res.result.type, profil_steps=res.result.profil_steps, dtype=np.float32)
		print("float32 drift: "+str(np.linalg.norm(np.array(res.result.end)-np.array(turn32.end))))
	return segments
		

paths, attributes = svgpathtools.svg2paths('test/track_layout.svg')

class TrackInfo():
	def __init__(self, name, category, width, author):
		self.name = name
		self.category = category
		self.author = author
fig = plt.figure()
#plt.gca().invert_yaxis()
ax = fig.add_subplot(111)

for svg_path in paths[1:2]:
	mat_path = convert_closed_path(svg_path)
	patch = patches.PathPatch(mat_path, facecolor='none', lw=1)
	ax.add_patch(patch)

for svg_path in paths[1:2]:
	mat_path = tangent_paths(svg_path)
	patch = patches.PathPatch(mat_path, facecolor='none', lw=1, edgecolor='red')
	ax.add_patch(patch)

start = timer()
segments = None
for svg_path in paths[1:2]:
	segments = fit_torcs_segments(svg_path, fig)
	#patch = patches.PathPatch(torcs_path, facecolor='none', lw=1, edgecolor='green')
	#ax.add_patch(patch)
end = timer()
print(end - start)

segments = fix_closing_for_float32(segments)

track_info = TrackInfo("Test1", "road", 8.0, "0x1af")
with open('test-1.xml', 'w') as f:
	f.write(torcs_trackxml.create_new(track_info, segments))

ax.autoscale(enable=True, axis='both', tight=None)
plt.axis('equal')
plt.show()
