#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""\
SVG.py - Construct/display SVG scenes.

The following code is a lightweight wrapper around SVG files. The metaphor
is to construct a scene, add objects to it, and then write it to a file
to display it.
"""

display_prog = 'display'  # Command to execute to display images.


class Scene:
	def __init__(self, name="svg", height=150, width=750):
		self.name = name
		self.items = []
		self.height = height
		self.width = width
		return

	def add(self, item):
		self.items.append(item)

	def strarray(self):
		var = ["<?xml version=\"1.0\"?>\n",
		       "<svg height=\"%d\" width=\"%d\" >\n" % (self.height, self.width),
		       " <g style=\"fill-opacity:1.0; stroke:black;\n",
		       "  stroke-width:1;\">\n"]
		for item in self.items: var += item.strarray()
		var += [" </g>\n</svg>\n"]
		return var

	def write_svg(self, filename=None):
		if filename:
			self.svgname = filename
		else:
			self.svgname = self.name + ".svg"
		file = open(self.svgname, 'w')
		file.writelines(self.strarray())
		file.close()
		return

	def display(self, prog=display_prog):
		os.system("%s %s" % (prog, self.svgname))
		return


class Line:
	def __init__(self, start, end):
		self.start = start  # xy tuple
		self.end = end  # xy tuple
		return

	def strarray(self):
		return ["  <line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" />\n" % \
		        (self.start[0], self.start[1], self.end[0], self.end[1])]


class Circle:
	def __init__(self, center, radius, color):
		self.center = center  # xy tuple
		self.radius = radius  # xy tuple
		self.color = color  # rgb tuple in range(0,256)
		return

	def strarray(self):
		return ["  <circle cx=\"%d\" cy=\"%d\" r=\"%d\"\n" % \
		        (self.center[0], self.center[1], self.radius),
		        "    style=\"fill:%s;\"  />\n" % colorstr(self.color)]


class TE:
	def __init__(self, origin, type=""):
		self.origin = [origin, 95]
		self.point_left = [(origin + 8), 80]
		self.point_right = [(origin - 8), 80]
		if type == "ins":
			self.fill = "#6677CC"
		elif type == "del":
			self.fill = "##DD88BB"
		else:
			self.fill = "none"
		return

	def strarray(self):
		return ["  <polygon points=\"%d,%d %d,%d %d,%d \" style=\"fill:%s\"/>\n" % \
		        (self.point_right[0], self.point_right[1], self.point_left[0], self.point_left[1], self.origin[0],
		         self.origin[1], self.fill)]


class Rectangle:
	def __init__(self, origin, height, width, color, opacity=None):
		self.origin = origin
		self.height = height
		self.width = width
		self.color = color
		if opacity != None:
			self.opacity = opacity
		else:
			self.opacity = 1
		return

	def strarray(self):
		return ["  <rect x=\"%d\" y=\"%d\" height=\"%d\"\n" % \
		        (self.origin[0], self.origin[1], self.height),
		        "    width=\"%d\" style=\"fill:%s;,fill-opacity:%s\" />\n" % \
		        (self.width, colorstr(self.color), self.opacity)]


class Text:
	def __init__(self, origin, text, size=6):
		self.origin = origin
		self.text = text
		self.size = size
		return

	def strarray(self):
		return ["  <text x=\"%d\" y=\"%d\" font-size=\"%d\" font-family=\"sans\" font-weight=\"100\" >\n" % \
		        (self.origin[0], self.origin[1], self.size),
		        "   %s\n" % self.text,
		        "  </text>\n"]


class Allele:
	def __init__(self, start, end, name=None, color_num=1):
		rectangle_width = end - start
		color_list = [[171, 138, 222], [221, 136, 187], [206, 146, 135], [222, 213, 138], [206, 146, 135]]
		self.rect = Rectangle([start, 95], 10, rectangle_width, color_list[color_num], 0.4)
		self.name = name
		self.start = start
		self.end = end

	def removeName(self):
		self.name = None
		return self

	def getWidth(self):
		return (self.end - self.start)

	def strarray(self):
		if self.name != None:
			self.label = Text([self.start, 117], self.name, 8)
			self.obj = self.rect.strarray() + self.label.strarray()
		else:
			self.obj = self.rect.strarray()
		return self.obj


def colorstr(rgb): return "#%x%x%x" % (int(rgb[0] / 16), int(rgb[1] / 16), int(rgb[2] / 16))
