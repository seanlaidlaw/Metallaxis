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
		# add 50 to width so that annotation doesn't get cut off
		var = ["<?xml version=\"1.0\"?>\n",
		       "<svg xmlns=\"http://www.w3.org/2000/svg\" height=\"%d\" width=\"%d\" >\n" % (self.height, (self.width + 50)),
                       "<style>",
                       ".te polygon {fill-opacity:0.7}",
                       ".te text {visibility: hidden; stroke-width: 0px; font-size: 8pt; font-family: sans-serif;}",
                       ".te:hover polygon {fill-opacity: 1;}",
                       ".te:hover text {visibility: visible;}",
                       ".allele rect {fill-opacity:0.7}",
                       ".allele text {visibility: hidden; stroke-width: 0px; font-size: 8pt; font-family: sans-serif;}",
                       ".allele:hover rect {fill-opacity: 1;}",
                       ".allele:hover text {visibility: visible;}",
                       "</style>",
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
	def __init__(self, origin, name, current_impact,current_annotation, type=""):
		self.origin = [origin, 95]
		self.name = name
		self.current_annotation = current_annotation
		self.point_right = [(origin + 8), 80]
		self.point_left = [(origin - 8), 80]
		if type == "ins":
			self.fill = "#6677CC"
		elif type == "del":
			self.fill = "#DD88BB"
		else:
			self.fill = "#fff"

		# set different colors & widths depending on annotation impact
		if current_impact == "HIGH":
			self.strokeWidth = "2"
			self.stroke = "red"
		elif current_impact == "MODERATE":
			self.strokeWidth = "2"
			self.stroke = "yellow"
		elif current_impact == "MODIFIER":
			self.strokeWidth = "1"
			self.stroke = "green"
		elif current_impact == "LOW":
			self.strokeWidth = "1"
			self.stroke = "green"
		else:
			self.strokeWidth = "1"
			self.stroke = "black"

		return

	def strarray(self):
		self.obj = ["  <polygon points=\"%d,%d %d,%d %d,%d \" style=\"fill:%s\" stroke=\"%s\" stroke-width=\"%s\" />\n" % (self.point_right[0], self.point_right[1], self.point_left[0], self.point_left[1], self.origin[0], self.origin[1], self.fill, self.stroke, self.strokeWidth)]

		if self.name != ".":
			self.label = Text([self.point_left[0], 70], self.name, 8).strarray()
		else:
			self.label = Text([self.point_left[0], 70], " ", 8).strarray()

		if self.current_annotation is not None:
			self.label += Text([self.point_left[0], 60], self.current_annotation, 8, fontColor="#2062ba", fontWeight=700).strarray()

		self.obj = ['<g class="te">'] + self.obj + self.label + ['</g>']
		return self.obj


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
	def __init__(self, origin, text, size=6, fontColor="#303030", fontStyle="normal", fontWeight=100):
		self.origin = origin
		self.text = text
		self.size = size
		self.fontColor = fontColor
		self.fontStyle = fontStyle
		self.fontWeight = fontWeight
		return

	def strarray(self):
		return ["  <text x=\"%d\" y=\"%d\" fill=\"%s\" font-style=\"%s\" font-size=\"%d\" font-family=\"sans\" font-weight=\"%d\" >\n" % \
		        (self.origin[0], self.origin[1], self.fontColor, self.fontStyle, self.size, self.fontWeight),
		        "   %s\n" % self.text,
		        "  </text>\n"]


class Allele:
	def __init__(self, start, end, name, biotype, description, color_num=1):
		rectangle_width = end - start
		color_list = [[171, 138, 222], [221, 136, 187], [206, 146, 135], [222, 213, 138], [206, 146, 135]]
		self.rect = Rectangle([start, 95], 10, rectangle_width, color_list[color_num], 0.4)
		self.name = name
		self.start = start
		self.end = end
		self.biotype = biotype
		self.description = description

	def getWidth(self):
		return (self.end - self.start)

	def strarray(self):
		self.label = Text([self.start, 117], self.name, 8, fontWeight=700, fontColor="#303030").strarray()
		if self.biotype != None:
		    self.label += Text([self.start, 127], self.biotype, 8, fontColor="#6c79c5").strarray()
		if self.description != None:
		    self.label += Text([self.start, 137], self.description, 8, fontColor="#787a88", fontStyle="italic").strarray()

		self.obj = ['<g class="allele">'] + self.rect.strarray() + self.label + ['</g>']
		return self.obj


def colorstr(rgb): return "#%x%x%x" % (int(rgb[0] / 16), int(rgb[1] / 16), int(rgb[2] / 16))
