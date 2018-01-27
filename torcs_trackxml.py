#! /usr/bin/python3
import codecs
from lxml import etree
from lxml.builder import E
import datetime
from pprint import pprint
import collections
import numpy as np
import pdb
from track_segment import Turn


#class Texture():
#	__init__(self, name, tex_type, size, mipmap):
#		pass
def filesys_friendly_str(name):
	return name.lower().replace(' ','-').strip()

def header_comment(filename, author, email):
	email = 'unknown'
	version = 'unknown'
	now_time = datetime.datetime.now()
	created = now_time.strftime("%I:%M%p on %B %d, %Y")
	copyright = "(C) %s by %s"%(now_time.strftime("%Y"), author)
	gen_cmd  = "TODO"
#	pprint(locals())
	lines = ('\nfile: {filename}\ncreated: {created}\n\ncopyright: {copyright}\n'
		'email: {email}\nversion: {version}\n initially generated with:\n  {gen_cmd} \n').format(**locals())
	return etree.Comment(lines)

def license_comment():
	license_lines = ('''\nThis program is free software; you can redistribute it and/or modify\n'''
		'''it under the terms of the GNU General Public License as published by\n'''
		'''the Free Software Foundation; either version 2 of the License, or\n'''
		'''(at your option) any later version.\n''')
	return etree.Comment(license_lines)

def SECTION(name, *children):
	return E.section({"name": name}, *children)

def ATTNUM(name, val, **opt):
	return E.attnum(collections.OrderedDict([ ('name', name), ('val', val) ] + list(opt.items()) ))

def ATTSTR(name, val, **opt):
	return E.attstr(collections.OrderedDict([ ('name', name), ('val', val) ] + list(opt.items()) ))

def track_header(name, category, trackformat_version, author, desc):
	return SECTION("Header", 
		ATTSTR("name", name),
		ATTSTR("category",category),
		ATTNUM("version", trackformat_version),
		ATTSTR("author", author),
		ATTSTR("description", desc)
	)
def surface(name, color1, color2, texture, friction, rolling_resistance, bump_tex, roughness):
	return SECTION(name) #TODO complete

def track_surfaces():
	return SECTION("Surfaces",
		etree.Entity("default-surfaces")
	)

def track_graphic(desc_file_3d):
	return SECTION("Graphic",
		ATTSTR("3d description", desc_file_3d),
		SECTION("Terrain Generation",
		      ATTNUM("track step", "20", unit="m"),
		      ATTNUM("border margin", "50", unit="m"),
		      ATTNUM("border step", "30", unit="m"),
		      ATTNUM("border height", "15", unit="m"),
		      ATTNUM("orientation", "counter-clockwise")
	))

def track_side(type_str, start_width, end_width, surface):
	return SECTION(type_str,
	      ATTNUM("start width", start_width, unit="m"),
	      ATTNUM("end width", end_width, unit="m"),
	      ATTSTR("surface", surface),		
	)

def track_boundary(type_str, width, height, surface, style):
	return SECTION(type_str,
	      ATTNUM("width", width, unit="m"),
	      ATTNUM("height", height, unit="m"),
	      ATTSTR("surface", surface),
	      ATTSTR("style", style)		
	)

def track_segment(name, segment):
	props = [
		ATTSTR("type", segment.type),
		ATTNUM("z start", "0.0", unit="m"),
		ATTNUM("z end", "0.0", unit="m")
	]
	if isinstance(segment, Turn):
		props.append(ATTNUM("radius", str(segment.radius), unit="m"))
		props.append(ATTNUM("end radius", str(segment.end_radius), unit="m"))
		props.append(ATTNUM("arc", str(segment.arc/np.pi*180.0), unit="deg"))
		props.append(ATTNUM("profil steps", str(int(segment.profil_steps))))
		banking_val = "5.0"
		if segment.type == 'lft':
			banking_val = "-5.0"
		#props.append(ATTNUM("banking start", banking_val, unit="deg"))
		#props.append(ATTNUM("banking end", banking_val, unit="deg"))

	props.append(track_side("Left Side", "5.0", "5.0", "grass"))
	props.append(track_side("Right Side", "5.0", "5.0", "grass"))
	return SECTION(name, *props)
	
def track_segments(segments):
	xml_segments = []
	for i, segment in enumerate(segments):
		xml_segments.append(track_segment(str(i), segment))

	return SECTION("Track Segments",
		*xml_segments
	)
	
	

def track_main(width, segments):
	return SECTION("Main Track",
		ATTNUM("width", str(width), unit="m"),
		ATTNUM("profil steps length", "4.0", unit="m"),
		ATTSTR("surface", "asphalt2-lines"),
		track_side("Left Side", "5.0", "5.0", "grass"),
		track_boundary("Left Border", "0.5", "0.05", "curb-5cm-r", "plan"),
		track_boundary("Left Barrier", "0.1", "1.0", "barrier", "curb"),
		track_side("Right Side", "5.0", "5.0", "grass"),
		track_boundary("Right Border", "0.5", "0.05", "curb-5cm-r", "plan"),
		track_boundary("Right Barrier", "0.1", "1.0", "barrier", "curb"),
		track_segments(segments)
	)

def create_new(track_info, segments):
	track_name = "Test 1"
	file_basename = filesys_friendly_str(track_name)
	category = "road"
	trackformat_version = "4"
	author = "Jana Doe"
	desc = "Test track"
	width = 9.0

	page = (
		E.params({"name": track_name, "type": "trackdef", "mode": "mw"},       # create an Element called "html"
			track_header(track_name, category, trackformat_version, author, desc ),
			track_surfaces(),
			track_graphic(file_basename+"-trk.ac"),
			track_main(width, segments)
		)
	)

	page.addprevious(header_comment(file_basename+".xml", author, "jane@doe.do"))
	page.addprevious(license_comment())

	tree = etree.ElementTree(page)
	doctype_str = ('<!DOCTYPE params SYSTEM "../../../../src/libs/tgf/params.dtd" [\n'
			'<!-- general definitions for tracks -->\n'
		'<!ENTITY default-surfaces SYSTEM "../../../data/tracks/surfaces.xml">\n'
		'<!ENTITY default-objects SYSTEM "../../../data/tracks/objects.xml">\n'
		']>')
	return etree.tostring(tree, pretty_print=True, xml_declaration=True, doctype=doctype_str, encoding='UTF-8').decode('utf-8')

def indent(elem, level=0):
    i = "\n" + level*"\t"
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + "\t"
        for e in elem:
            indent(e, level+1)
            if not e.tail or not e.tail.strip():
                e.tail = i + "\t"
        if not e.tail or not e.tail.strip():
            e.tail = i
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = i

def update(file_name, segments):
	reader = codecs.open(file_name, 'r',
                     encoding='utf-8',
                     errors='replace')
	lines = reader.readlines()

	parser = etree.XMLParser(ns_clean=True, resolve_entities=False, remove_blank_text=False)
	tree = etree.fromstring(''.join(lines).encode('utf8', 'replace'), parser)#.getroottree()
	preceding_elem = []
	'''for elem in tree.getroot().itersiblings(preceding=True):
		print(elem.text)
		print(elem.tail)
		preceding_elem.append(elem)
	for elem in reversed(preceding_elem):
		elem.tail = 'Wassup'
		tree.getroot().addprevious(elem)
	for elem in tree.getroot().itersiblings(preceding=True):
		print(elem.text)
		print(elem.tail)
		preceding_elem.append(elem)
	pdb.set_trace()
	'''
	#for child in tree.getroot():
	#	print('tag: {} sourceline:{}'.format(child.tag, child.sourceline))
	#print(tree.getroot().sourceline)
	'''for elem in tree.getroot().iter():
		if elem.tail:
			elem.tail = elem.tail.replace("\n", "", 1)
		if not type(elem) is etree._Entity and elem.text:
			str_text = elem.text
			print(str_text)
			elem.text = str_text.replace("\n", "", 1)
		
		pdb.set_trace()'''
	#tree = etree.fromstring(lines.encode('utf8', 'replace')).getroottree()
	#print(tree)
	#print(etree.tostring(tree.getroottree()))
	trackSection = tree.getroottree().find('section[@name="Main Track"]')
	trackWidth = float(trackSection.find('attnum[@name="width"]').get('val'))
	defaultStepLength = float(trackSection.find('attnum[@name="profil steps length"]').get('val'))
	print(trackWidth)
	tree.addprevious(license_comment())
	segSection = trackSection.find('section[@name="Track Segments"]')
	segs = segSection.getchildren()
	for seg in segs:
		if seg.tag is etree.Comment:
			continue
		seg_name = seg.get('name')
		#pdb.set_trace()
		seg.set('name', seg_name+"olaf was here")

	with open("mod_"+file_name, 'w') as f:
		_str = etree.tostring(tree, xml_declaration=False, pretty_print=False, encoding='UTF-8').decode('utf-8')
		#_str = _str.replace('\n\n', '\n') # try to preserve the original format, to_string auto-adds '\n' newlines, remove them again
		headerend_line = tree.sourceline-1
		f.write(''.join(lines[:headerend_line]))
		f.write(_str)
		pdb.set_trace()
		f.write(tree.tail)
		treeend_line = _str.count('\n') + headerend_line+1
		f.write(''.join(lines[treeend_line:]))
			
	return

if __name__ == "__main__":
	# execute only if run as a script
	
	update("alpine-1.xml", None)

