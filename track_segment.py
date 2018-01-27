class Turn():
	def __init__(self, _type, start, end, arc, verts, profil_steps):
		self.type = _type
		self.start = start
		self.end = end
		self.arc = arc
		self.verts = verts
		self.profil_steps = profil_steps

	def set_radius(self, start_radius, end_radius):
		self.radius = start_radius
		self.end_radius = end_radius
