#!/usr/bin/env python

#--------------------------------------------------------------
# Function to add the aperture class instances to the linac lattice,
# and a function to create the array with aperture nodes and losses.
#--------------------------------------------------------------

import math
import sys
import os

from orbit.py_linac.lattice import LinacApertureNode
from orbit.py_linac.lattice import Quad
from orbit.py_linac.overlapping_fields import OverlappingQuadsNode

def Add_quad_apertures_to_lattice(accLattice):
	"""
	Function will add Aperture nodes at the entrance and exit of quads.
	It returns the list of Aperture nodes.
	"""
	aprtNodes = []
	node_pos_dict = accLattice.getNodePositionsDict()
	quads = accLattice.getNodesOfClasses([Quad,OverlappingQuadsNode])
	for node in quads:
		if(isinstance(node,Quad)):
			if(node.hasParam("aperture") and node.hasParam("aprt_type")):
				shape = node.getParam("aprt_type")
				a = node.getParam("aperture")
				node_name = node.getName()
				(posBefore, posAfter) = node_pos_dict[node]
				apertureNodeBefore = LinacApertureNode(shape,a/2.0,a/2.0,posBefore)
				apertureNodeAfter = LinacApertureNode(shape,a/2.0,a/2.0,posAfter)
				apertureNodeBefore.setName(node_name+":AprtIn")
				apertureNodeAfter.setName(node_name+":AprtOut")
				apertureNodeBefore.setSequence(node.getSequence())
				apertureNodeAfter.setSequence(node.getSequence())
				node.addChildNode(apertureNodeBefore,node.ENTRANCE)
				node.addChildNode(apertureNodeAfter,node.EXIT)
				aprtNodes.append(apertureNodeBefore)
				aprtNodes.append(apertureNodeAfter)
		if(isinstance(node,OverlappingQuadsNode)):
			# place to add aperture for the overlapping fields quads
			nParts = node.getnParts()
			simple_quads = node.getQuads()
			quad_centers = node.getCentersOfField()
			node_start_pos = node_pos_dict[node][0]
			pos_part_arr = []
			s = 0.
			for part_ind in range(nParts):
				pos_part_arr.append(s)
				s += node.getLength(part_ind)
			for quad_ind in range(len(simple_quads)):
				quad = simple_quads[quad_ind]
				shape = quad.getParam("aprt_type")
				a = quad.getParam("aperture")
				quad_name = quad.getName()				
				length = quad.getLength()
				pos_z = quad_centers[quad_ind]
				posBefore = pos_z - length/2.
				posAfter = pos_z + length/2.
				for part_ind in range(nParts-1):
					if(posBefore >= pos_part_arr[part_ind] and posBefore < pos_part_arr[part_ind+1]):
						apertureNodeBefore = LinacApertureNode(shape,a/2.0,a/2.0,posBefore + node_start_pos)
						apertureNodeBefore.setName(quad_name+":AprtIn")
						apertureNodeBefore.setSequence(node.getSequence())
						node.addChildNode(apertureNodeBefore,node.BODY,part_ind,node.BEFORE)
						aprtNodes.append(apertureNodeBefore)
					if(posAfter > pos_part_arr[part_ind] and posAfter <= pos_part_arr[part_ind+1]):
						apertureNodeAfter = LinacApertureNode(shape,a/2.0,a/2.0,posAfter + node_start_pos)
						apertureNodeAfter.setName(quad_name+":AprtOut")
						apertureNodeAfter.setSequence(node.getSequence())
						node.addChildNode(apertureNodeAfter,node.BODY,part_ind,node.AFTER)
						aprtNodes.append(apertureNodeAfter)
	return aprtNodes
	
def GetLostDistributionArr(aprtNodes,bunch_lost):
	"""
	Function returns the array with [aptrNode,sum_of_losses]
	The sum_of_losses is a number of particles or the sum of macrosizes if the 
	particle attribute "macrosize" is defined.
	"""
	lossDist_arr = []
	aprtPos_arr = []
	#--- first we will sort apertures according to the position
	aprtNodes = sorted(aprtNodes, key = lambda x: x.getPosition(), reverse = False)
	for aprt_node in aprtNodes:
		aprtPos_arr.append(aprt_node.getPosition())
		loss_sum = 0.
		lossDist_arr.append([aprt_node,loss_sum])
	if(len(aprtPos_arr) <= 0): return lossDist_arr
	#-----------------------------------------------------------------
	def indexFindF(pos_arr,pos,ind_start = -1,ind_stop = -1):
		""" This function will find the index of nearest to pos point in  pos_arr"""
		if(ind_start < 0 or ind_stop < 0): 
			ind_start = 0 
			ind_stop = len(pos_arr) - 1
			return indexFindF(pos_arr,pos,ind_start,ind_stop)
		if(abs(ind_start - ind_stop) <= 1):
			dist0 = abs(pos_arr[ind_start] - pos)
			dist1 = abs(pos_arr[ind_stop] - pos)			
			if(dist0 <= dist1): return ind_start
			return ind_stop
		ind_mdl = int((ind_start + ind_stop)/2.0)
		if(pos_arr[ind_start] <= pos and pos <pos_arr[ind_mdl]):
			return indexFindF(pos_arr,pos,ind_start,ind_mdl)
		else:
			return indexFindF(pos_arr,pos,ind_mdl,ind_stop)
	#-----------------------------------------------------------------
	has_macrosize_partAttr = bunch_lost.hasPartAttr("macrosize")
	if(not bunch_lost.hasPartAttr("LostParticleAttributes")): return lossDist_arr
	macroSize = bunch_lost.macroSize()
	nParticles = bunch_lost.getSize()
	for ind in range(nParticles):
		pos = bunch_lost.partAttrValue("LostParticleAttributes",ind,0)
		if(pos < aprtPos_arr[0] or pos > aprtPos_arr[len(aprtPos_arr)-1]): continue
		pos_ind = indexFindF(aprtPos_arr,pos)
		if(has_macrosize_partAttr):
			macroSize = bunch_lost.partAttrValue("macrosize",ind,0)
			lossDist_arr[pos_ind][1] += macroSize
			continue
		lossDist_arr[pos_ind][1] += 1.0
	return lossDist_arr
			

def AddScrapersAperturesToLattice(accLattice,node_name,x_size,y_size,aprtNodes):
	"""
	Function will add the rectangular Aperture node (shape=3) at the node with a particular name.
	Parameters x_size and y_size are full horizontal and vertical sizes of the aperture.
	"""
	shape = 3
	node_pos_dict = accLattice.getNodePositionsDict()
	node = accLattice.getNodesForName(node_name)[0]
	(posBefore, posAfter) = node_pos_dict[node]
	apertureNode = LinacApertureNode(shape,x_size/2.0,y_size/2.0,posBefore)
	apertureNode.setName(node_name+":Aprt")
	apertureNode.setSequence(node.getSequence())
	node.addChildNode(apertureNode,node.ENTRANCE)
	aprtNodes.append(apertureNode)
	aprtNodes = sorted(aprtNodes, key = lambda x: x.getPosition(), reverse = False)
	return aprtNodes			
		
			
		
		
	