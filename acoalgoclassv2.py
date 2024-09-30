from qgis.core import QgsVectorLayer, QgsApplication, QgsProject, QgsField, QgsFeature, QgsGeometry, QgsPointXY, QgsPoint, edit, QgsWkbTypes
from qgis.core import QgsSpatialIndex, QgsVectorFileWriter
from qgis.gui import QgsMapCanvas, QgsMapTool, QgsLayerTreeMapCanvasBridge
from PyQt5.QtWidgets import QApplication
from PyQt5.QtGui import QColor
from qgis.PyQt.QtCore import QTimer, Qt
from PyQt5.QtCore import QVariant, QThread
from collections import OrderedDict
import queue
import numpy as np
import shutil, os
import sys

class AntColonyOptimization():
        """Ant Colony Optimization is an algorithm imitating the behavior of ants in terms of reaching point to point
        \nAnts = Number of Ants
        \nQ = Rate of deposit of pheromones on each iteration
        \nAlpha = Value importance of Pheromone Values
        \nBeta = Value importance in terms of distance and other factors
        \nroadLayer = Vector layer of your road
        \nsiteNodes = Vector layer of sites need to be visited
        \nroadLinks = Vector layer of route links for intersections
        """
        #CODE FUNCTION PARTIAL COMPLETE
        def __init__(self, Ants:int, Q:int, alpha:float, beta:float, rho:int, roadLayer:QgsVectorLayer, siteNodes:QgsVectorLayer, roadObstacles:QgsVectorLayer):
            print("INIT ACO")
            self.ants = Ants
            self.Q = Q
            self.alpha = alpha
            self.beta = beta
            self.rho = rho
            self.pheromoneLayer = None
            self.roadLayer = roadLayer
            #self.roadDirections = None
            self.siteNodes = siteNodes
            self.roadObstacles = roadObstacles
            #constant values OSM BASED VALUES
            self.roadType=("primary", "secondary", "tertiary", "residential", "unclassified")
            self.roadLinks=("primary_link","secondary_link","tertiary_link","lane_link","generalroad_link")
            self.osmType=("way")
            self.osm_Access = ("private","permissive", "no")
            #number of It
            self.num_iterations = 1
            #Simulation Point
            self.simullayer = QgsVectorLayer('Point?crs=EPSG:3857', 'Current Position', 'memory')
            self.providersim = self.simullayer.dataProvider()
            self.providersim.addAttributes([QgsField('id', QVariant.Int)])
            self.simullayer.updateFields()
            self.currpoint = None
            #RUNNING TIME
            self.iter_tour = 0
            self._stop_flag = False
            
        def initLayer(self):
            QgsProject.instance().addMapLayer(self.simullayer)
            
        def update_position(self):
            self.simullayer.dataProvider().truncate()
            
            # Create a new feature
            feature = QgsFeature()
            feature.setGeometry(QgsGeometry.fromPointXY(self.currpoint))
            feature.setAttributes([1])  # Set an id attribute
            self.providersim.addFeature(feature)
            
            # Update the layer
            self.simullayer.updateExtents()
            self.canvas.refresh()
            
        #INTERSECTION FUNCTION
        def pointIntersect(self, point:QgsPointXY, id:int, layer:QgsVectorLayer, flags = 1):
            spatial_index = QgsSpatialIndex(layer.getFeatures())
            potential_ids = spatial_index.nearestNeighbor(point, 10, 0.5)
            intersected_ids = []
            point_geo = QgsGeometry.fromPointXY(point)
            point_geo_buffer = point_geo.buffer(1e-3, 5)
            # Count the actual intersections
            for feature_id in potential_ids:
                if feature_id == id:
                    continue
                feature = layer.getFeature(feature_id)
                if feature.attribute("highway") in self.roadType or feature.attribute("highway") in self.roadLinks:
                    if feature.geometry().intersects(point_geo_buffer):
                        intersected_ids.append(feature_id)
            if flags == 1:
                return len(intersected_ids)
            elif flags == 2:
                return intersected_ids
            
        def returnPheromoneEdgeID(self, route:list, flags = 1, numneighbors = 10):
            spatial_index = QgsSpatialIndex(self.pheromoneLayer.getFeatures())
            geo, geobuff = None, None
            if flags == 1:
                geo = QgsGeometry.fromPolylineXY(QgsPointXY(p) for p in route)
                geobuff = geo.buffer(1e-1, 5)
            else:
                geo = QgsGeometry.fromPointXY(QgsPointXY(route[0]))
                geobuff = geo.buffer(1e-1, 5)
            potential_ids = spatial_index.nearestNeighbor(geo, numneighbors, 1)
            intersected_ids = []
            for feature_id in potential_ids:
                feature = self.pheromoneLayer.getFeature(feature_id)
                if flags == 1:
                    if feature.geometry().within(geobuff):
                        intersected_ids.append(feature_id)
            return intersected_ids
        
        #PHEROMONE FUNCTIONS    
        #Function to get edges on each nodes, going to be used for calculating pheromones
        def InitPheromones(self):
            """Creates A Temporary Layer Using your Road Layer as a reference to initialize edges with specified pheromone values
            """
            if self.roadLayer is None:
                print("You cannot Pass a None Layer!")
                return False
            if self.loadpheromoneLayer():
                print('Pheromone Data was Detected')
            else:
                print("Init Pheromones")
                self.pheromoneLayer = QgsVectorLayer("LineString?crs=epsg:3857", "pheromones", "memory")
                provider = self.pheromoneLayer.dataProvider()
                provider.addAttributes([QgsField("pheromone", QVariant.Double)])
                provider.addAttributes([QgsField("bpheromone", QVariant.Double)])
                provider.addAttributes([QgsField("distance", QVariant.Double)])
                provider.addAttributes([QgsField("feid", QVariant.Int)])
                road_type_field = QgsField("road_type", QVariant.String, len=254)
                provider.addAttributes([road_type_field])
                road_type_field = QgsField("status", QVariant.String, len=254)
                provider.addAttributes([road_type_field])
                self.pheromoneLayer.updateFields()
                
                #Road Layer and intersection directions
                with edit(self.pheromoneLayer):
                    print("CREATE PHEROMONE LINKS FOR ROAD LAYER")
                    for feature in self.roadLayer.getFeatures():
                        if feature.attribute("highway") not in self.roadType:
                            continue
                        geom = feature.geometry()
                        if geom.isMultipart():
                            for line in geom.asMultiPolyline():
                                for vertex in range(1, len(line)):
                                    point = []
                                    sv = line[vertex - 1]   #starting vertex
                                    ev = line[vertex]       #end vertex
                                    point.append(sv)        #add them to a list
                                    point.append(ev)
                                    edge = QgsGeometry.fromPolylineXY(QgsPointXY(p) for p in point)
                                    distance = edge.length()
                                    edgef = QgsFeature(self.pheromoneLayer.fields())
                                    edgef.setGeometry(edge)                                         
                                    edgef.setAttributes([
                                            1.0,
                                            1.0,
                                            distance,
                                            feature.id(),
                                            feature.attribute("highway"),
                                            "open"
                                    ])                                     
                                    self.pheromoneLayer.addFeature(edgef)
                        else:
                            for line in geom.asPolyline():
                                for vertex in range(1, len(line)):
                                    point = []
                                    sv = line[vertex - 1]   #starting vertex
                                    ev = line[vertex]       #end vertex
                                    point.append(sv)        #add them to a list
                                    point.append(ev)
                                    edge = QgsGeometry.fromPolylineXY(QgsPointXY(p) for p in point)
                                    distance = edge.length()
                                    edgef = QgsFeature(self.pheromoneLayer.fields())        
                                    edgef.setGeometry(edge)                                         
                                    edgef.setAttributes([
                                            1.0,
                                            1.0,
                                            distance,
                                            feature.id(),
                                            feature.attribute("highway"),
                                            "open"
                                    ])                                     
                                    self.pheromoneLayer.addFeature(edgef)
                    print("DONE!")
                    print("CREATE PHEROMONE LINKS FOR ROAD DIRECTIONS")
                    for feature in self.roadLayer.getFeatures():
                        if feature.attribute("highway") not in self.roadLinks:
                            continue
                        geom = feature.geometry()
                        if geom.isMultipart():
                            for line in geom.asMultiPolyline():
                                points = []  # Collect all points for the edge
                                for vertex in line:
                                    points.append(vertex)
                                if points:
                                    edge = QgsGeometry.fromPolylineXY(points)
                                    distance = edge.length()
                                    edgef = QgsFeature(self.pheromoneLayer.fields())
                                    edgef.setGeometry(edge)
                                    edgef.setAttributes([
                                            1.0,
                                            1.0,
                                            distance,
                                            feature.id(),
                                            feature.attribute("highway"),
                                            None
                                        ])
                                    self.pheromoneLayer.addFeature(edgef)
                        else:
                            # Handle single part geometries if needed
                            line = geom.asPolyline()
                            points = []  # Collect all points for the edge
                            for vertex in line:
                                points.append(vertex)
                            if points:  # Ensure there are points collected
                                edge = QgsGeometry.fromPolylineXY(points)
                                distance = edge.length()
                                edgef = QgsFeature(self.pheromoneLayer.fields())
                                edgef.setGeometry(edge)
                                edgef.setAttributes([
                                            1.0,
                                            1.0,
                                            distance,
                                            feature.id(),
                                            feature.attribute("highway"),
                                            None
                                ])
                                self.pheromoneLayer.addFeature(edgef)
                    print("DONE!")
            return True
        def returnPheromoneLayer(self):
            return self.pheromoneLayer
        def deletePheromoneLayer(self):
            directory = "data/pherolayer"
            # Check if the directory exists before attempting to delete it
            if os.path.exists(directory):
                shutil.rmtree(directory)
                print(f"{directory} and all its contents have been deleted")
            else:
                print(f"{directory} does not exist")
                return False
            return True
        def savepheromoneLayer(self):
            directory_path = "data/pherolayer"
            os.makedirs(directory_path, exist_ok=True)
            save_path = "data/pherolayer/pherolayer.shp"
            error = QgsVectorFileWriter.writeAsVectorFormat(self.pheromoneLayer, save_path, "UTF-8", self.pheromoneLayer.crs(), "ESRI Shapefile")
            if os.path.exists(save_path):
                print(f"Pheromone Layer already exists")
                return 1
            if error[0] == QgsVectorFileWriter.NoError:
                print("Pheromone Data was saved")
                return error[0]
            else:
                print(f"Error Saving Pheromone Data {error}")
                return error[0]
        def loadpheromoneLayer(self):
            save_path = "data/pherolayer/pherolayer.shp"
            self.pheromoneLayer = QgsVectorLayer(save_path, "pheromones", "ogr")
            if not self.pheromoneLayer.isValid():
                print("Layer was not loaded")
                return False
            else:
                print("Layer was sucessfully loaded")
                return True
        def resetpheroValues(self):
            if self.pheromoneLayer.isValid():
                provider = self.pheromoneLayer.dataProvider()
                attributeindex1 = self.pheromoneLayer.fields().indexFromName("pheromone")
                attributeindex2 = self.pheromoneLayer.fields().indexFromName("bpheromone")
                nvalue = 1.0
                return True
            else:
                print("Something wrong with your Pheromone Layer...")
                return False
        #UPDATE PHEROMONES
        def updatePheromones(self, id:list):
            """Function to Update Pheromone Layer
            \nArgs:
            \n  id = List of Feature Id's to be updated
            """
            with edit(self.pheromoneLayer):
                for feature in self.pheromoneLayer.getFeatures():
                    if feature.id() in id:
                        geom = feature.geometry()
                        pheromone_value = feature['pheromone']
                        pheromone_value *= (1 - self.rho)  # Evaporation
                        pheromone_value += self.Q / geom.length()  # Deposit
                        self.pheromoneLayer.changeAttributeValue(feature.id(),0, pheromone_value)
                    else:
                        pheromone_value = feature['pheromone']
                        pheromone_value *= (1 - self.rho)  # Evaporation
                        self.pheromoneLayer.changeAttributeValue(feature.id(),0, pheromone_value)
            self.pheromoneLayer.commitChanges()
            self.pheromoneLayer.triggerRepaint()
        #INTERNAL FUNCTIONS
        def detectNearbyRoadObstacles(self, id:int, point:QgsPointXY):
            """Detects Nearby Road Obstacles and see if it intersects

            Args:
                point (QgsGeometry): Current Point of where the next point will be to determine road obstacles
            Returns:
                bool: True if point is accessable, False if the point is not
            """
            spatial_index = QgsSpatialIndex(self.roadObstacles.getFeatures())
            potential_ids = spatial_index.nearestNeighbor(point, 10, 5)
            point_geo = QgsGeometry.fromPointXY(point)
            point_geo_buffer = point_geo.buffer(1, 5)
            #print("OBSTALCES: ",potential_ids)
            # Count the actual intersections
            for feature_id in potential_ids:
                if feature_id == id:
                    continue
                feature = self.roadObstacles.getFeature(feature_id)
                if feature.geometry().intersects(point_geo_buffer):
                    #print("YES")
                    if not feature.attribute("barrier") is None:
                        #print("BARRIER")
                        return False
                    if feature.attribute("access") in self.osm_Access:
                        #print("ACCESS")
                        return False
                    if not feature.attribute("noexit") is None:
                        #print("NOEXIT")
                        return False
            #print("reached obstacle searcher")
            return True
        def detectNearbySiteVisit(self, id:int, point:QgsPointXY, sitelist:list, straysiteid=[], operate=False):
            spatial_index = QgsSpatialIndex(self.siteNodes.getFeatures())
            potential_ids = spatial_index.nearestNeighbor(point, 10, 5)
            point_geo = QgsGeometry.fromPointXY(point)
            point_geo_buffer = point_geo.buffer(1e-2, 5)
            i=0
            for feature_id in potential_ids:
                feature = self.siteNodes.getFeature(feature_id)
                if feature_id == id:
                    print(feature_id)
                    if feature.geometry().intersects(point_geo_buffer) and feature_id in sitelist:
                        print("SITE DETECTED")
                        if operate:
                            print(f"{feature_id} {sitelist}")
                            sitelist.remove(feature_id)
                            print(f"{sitelist}")
                        return True
                else:
                    if feature.geometry().intersects(point_geo_buffer) and feature_id in sitelist:
                        print("STRAY SITE DETECTED")
                        if operate:
                            print(f"{feature_id} {sitelist}")
                            siteindex = sitelist.index(feature_id)
                            straysiteid.append({
                                "siteid":feature_id,
                                "rank":siteindex
                            })
                            sitelist.remove(feature_id)
                            print(f"{sitelist}")
                        return True
            return False
        def returnHeuristicInfo(self, id:list, point:QgsPointXY, sitelist, straysiteid):
            distance = 0.0
            for ids in id:
                feature = self.pheromoneLayer.getFeature(ids)
                distance += feature.attribute("distance")
                if self.detectNearbySiteVisit(None,point,sitelist,straysiteid):
                    distance = distance * 100
            return distance
        
        def continousTraversal(self,feature:QgsFeature,point:QgsPointXY,route:list,featureids:set,targetFeatureID:int,sitelist:list,progress=None,canvasup=None):
            line_geo = feature.geometry()
            line_part = None
            if line_geo.isMultipart():
                line_part = line_geo.asMultiPolyline()
            else:
                line_part = line_geo.asPolyline()
            previous_point = point
            current_point = point
            infoforclosestvertex = line_geo.closestVertex(previous_point)
            
            #print("ROUTE", route)
            if line_geo is None:
                print("Empty Feature Geometry")
                return None
            if infoforclosestvertex[3] != -1 and line_part[0][infoforclosestvertex[3]] not in route:
                #print("right")
                introute = []
                introute.append(previous_point)
                _skip = False
                for i in range(infoforclosestvertex[1]+1, len(line_part[0])):
                    current_point = line_part[0][i]
                    #print("CRPOINT RIGHT: ",current_point, "CRPOINT PREV: ",previous_point)
                    if current_point not in route and self.detectNearbyRoadObstacles(None, current_point):
                        #print("unique node")
                        introute.append(current_point)
                        self.iter_tour += 1
                    else:
                        #print("not unique")
                        current_point = line_part[0][i-1]
                        break
                    if self.detectNearbySiteVisit(targetFeatureID, current_point,sitelist):
                        #print("site visited")
                        break
                    if self.pointIntersect(current_point, feature.id(), self.roadLayer):
                        #print("intersection")
                        break
                #print("time to extend route")
                if len(introute) == 1:
                    return None
                else:
                    route.extend(introute)
                #print("introute",introute)
                fetid = self.returnPheromoneEdgeID(introute, numneighbors=100)
                #print("fids",fetid)
                featureids.extend(fetid)
                return current_point
            
            elif infoforclosestvertex[2] != -1 and line_part[0][infoforclosestvertex[2]] not in route:
                #print("left")
                introute = []
                introute.append(previous_point)
                _skip = False
                for i in range(infoforclosestvertex[1]-1, -1, -1):
                    current_point = line_part[0][i]
                    #print("CRPOINT LEFT: ",current_point, "CRPOINT PREV: ",previous_point)
                    if current_point not in route and self.detectNearbyRoadObstacles(None, current_point):
                        #print("unique node")
                        introute.append(current_point)
                        self.iter_tour += 1
                    else:
                        #print("not unique")
                        current_point = line_part[0][i+1]
                        break
                    if self.detectNearbySiteVisit(targetFeatureID, current_point,sitelist):
                        #print("site visited")
                        break
                    if self.pointIntersect(current_point, feature.id(), self.roadLayer)>1:
                        #print("intersection")
                        break
                #print("time to extend route")
                if len(introute) == 1:
                    return None
                else:
                    route.extend(introute)
                #print("introute",introute)
                fetid = self.returnPheromoneEdgeID(introute, numneighbors=100)
                #print("fids",fetid)
                featureids.extend(fetid)
                return current_point
            
            print("???")
            
            return None
        def findNodeChoices(self, feature:QgsFeature, point:QgsPointXY, route:list, targetFeatureID:int, directionTolerance, backtracking=False):
            """Finds available node choices

            Args:
                feature (QgsFeature object): Feature line for to traverse through
                point (QgsGeometry:QgsPointXY object): Point as a reference
                route (list): List of QgsPointXY objects which are visited
                feature_ids: ID of edges to extract from pheromone layer, will use this as a reference to calculate pheromone values
            
            Returns:
            list = list:QgsPointXY = list of QgsPointXY
            """
            line_geo = feature.geometry()
            line_part = None
            if line_geo.isMultipart():
                line_part = line_geo.asMultiPolyline()
            else:
                line_part = line_geo.asPolyline()
            current_point = point
            infoforclosestvertex = line_geo.closestVertex(current_point)
            
            #print("NODE CHOICES\n\n")
            #[0] = QgsPointObject of the closest node [1] = index of the node [2] previous node [3] next node [4] distance of point
            # Check if the closest vertex has a valid previous or next node
            if not feature.attribute("osm_type") in self.osmType:
                return None
            if feature.attribute("oneway") == 'yes' and feature.attribute("highway") in self.roadType:
                #print("oneway")
                if feature.attribute("access") not in self.osm_Access:
                    if infoforclosestvertex[3] != -1:
                        #print("yes")
                        pointt = line_part[0][infoforclosestvertex[1]+1]
                        if pointt not in route and self.detectNearbyRoadObstacles(None, pointt):
                            return pointt   #QGSPOINT RETURN
                return None

            elif feature.attribute("osm_type") in self.osmType and feature.attribute("highway") in self.roadLinks:
                #rint("road_link")
                introute=[]
                introute.append(current_point)
                if infoforclosestvertex[3] != -1:
                    #print("yes")
                    for i in range(infoforclosestvertex[1]+1, len(line_part[0])):
                        point1 = line_part[0][i]
                        if not self.detectNearbyRoadObstacles(None, point1) or point1 in route:
                            return None
                        introute.append(point1)
                    if len(introute)==len(line_part[0]):
                        introute1 = []
                        introute1.append(introute)
                        directionTolerance[0] = True   #this determines that the route given was a road link
                        return introute1               #LIST OF QGS POINTS
                return None
            else:
                #print("two way")
                if feature.attribute("osm_type") in self.osmType and feature.attribute("access") not in self.osm_Access and feature.attribute("highway") in self.roadType:
                    nodes = []
                    if infoforclosestvertex[3] != -1:
                        #print("right")
                        if line_part[0][infoforclosestvertex[3]] not in route and self.detectNearbyRoadObstacles(None, line_part[0][infoforclosestvertex[3]]):
                            #print("yes")
                            nodes.append(line_part[0][infoforclosestvertex[3]])
                    if infoforclosestvertex[2] != -1:
                        #print("left")
                        if line_part[0][infoforclosestvertex[2]] not in route and self.detectNearbyRoadObstacles(None, line_part[0][infoforclosestvertex[2]]):
                            #print("yes")
                            nodes.append(line_part[0][infoforclosestvertex[2]])               
                    return nodes                                        #LIST OF QGS POINTS
        
        # Define a function to run the ACO algorithm
        #ant traversal function
        def ant_traverse(self, starting_node:QgsPointXY, ant:int, site_tovisit:list, starting_id, progress_callback=None, ):
            print(f"ANT {ant}... Traversing....")
            current_node = starting_node    #point xy object
            target_site = site_tovisit[0]   #id
            _CONST_NUMBERSITES = len(site_tovisit)   #number of sites to visit before backtracking
            num_sites_visited = []
            _CONST_SITES = tuple(site_tovisit.copy())
            iterate = 0     #target site rank
            
            #FEATURE ID FOR EACH VISITED ROUTE
            feature_id_for_pheromone = []
            #Route
            routesites = []
            route = []
            #list of routes, all sites and backtracking
            routes= []
            pheromoneIDS = []
            #InitialDecision
            InitialDecision = True
            directionTol = [False]
            backtracking = False
            backpoint = None
            self.iter_tour = 0
            
            #for stray sites
            stray_site_id = []
            
            #HARD RESET FOR ANTS WHO GOT STUCK
            HARD_RESET = 1
            HARD_RESET_MAX = 100
            
            while True:
                #print("ROUTE: ",route)
                if self._stop_flag:
                    return route, None
                if current_node is None:
                    #print("NODE IS NONE")
                    HARD_RESET += HARD_RESET * 1.15
                    #print(HARD_RESET)
                    if len(routesites)>0 and HARD_RESET > HARD_RESET_MAX and not backtracking:
                        print("\n\n\nHARD RESET with checkpoints")
                        temp = routesites.pop()
                        #print(temp)
                        num_sites_visited.remove(temp['SITEID'])
                        site_tovisit.insert(temp["INDEX"],temp['SITEID'])
                        #print(f"SITETOVISIT:{num_sites_visited}")
                        
                        if temp['INDEX'] < iterate:
                            target_site = _CONST_SITES[temp['INDEX']]
                            iterate = temp['INDEX']
                        
                        if routesites:
                            temp = routesites[-1]
                            route = temp["ROUTE"].copy()
                            feature_id_for_pheromone = temp["PHEROIDS"].copy()
                            current_node = route[-1]
                        else:
                            print("big hard reset")
                            route, routesites = [], []
                            feature_id_for_pheromone = []
                            num_sites_visited=[]
                            site_tovisit = list(_CONST_SITES)
                            current_node = starting_node
                            target_site = _CONST_SITES[0]
                        HARD_RESET = 1
                        continue
                    elif len(routesites)>0 and not backtracking:
                        print("reset")
                        temp = routesites[-1]
                        route = temp["ROUTE"].copy()
                        feature_id_for_pheromone = temp["PHEROIDS"].copy()
                        current_node = route[-1]
                        continue
                    else:
                        print("big hard reset")
                        route, routesites = [], []
                        feature_id_for_pheromone = []
                        num_sites_visited=[]
                        site_tovisit = list(_CONST_SITES)
                        current_node = starting_node
                        target_site = _CONST_SITES[0]
                        iterate = 0
                        HARD_RESET = 1
                        continue
                        
                    
                print(f"Current Node:{current_node} ITER TOUR: {self.iter_tour}")
                if progress_callback:
                    progress_callback.emit({
                        'route':route
                    })
                    
                route = list(OrderedDict.fromkeys(route))
                if self.detectNearbySiteVisit(target_site, current_node, site_tovisit, stray_site_id, operate = True):
                    #print("DETECT STORE ROUTES")
                    print(target_site)
                    HARD_RESET = 1
                    if stray_site_id:
                        tempid = stray_site_id.pop()
                        num_sites_visited.append(tempid['siteid'])
                        routesites.append({"ROUTE":route.copy(), 
                                        "PHEROIDS":feature_id_for_pheromone.copy(),
                                        "SITEID":tempid['siteid'],
                                        "INDEX":tempid['rank']})
                        #rint("ROUTESITES",routesites)
                    else:
                        num_sites_visited.append(target_site)
                        routesites.append({"ROUTE":route.copy(), 
                                        "PHEROIDS":feature_id_for_pheromone.copy(),
                                        "SITEID":target_site,
                                        "INDEX":iterate})
                        iterate += 1
                        #print("ROUTESITES",routesites)
                        if site_tovisit:
                            target_site = _CONST_SITES[iterate]
                        
                    if len(num_sites_visited) == _CONST_NUMBERSITES and not backtracking:
                        #print("DONE, BACKTRACKING")
                        HARD_RESET = 0.1
                        #print("EST",routesites)
                        backtracking = True
                        backpoint = route[len(route)-2]
                        current_node = route[-1]
                        
                        p1 = routesites[-1]['ROUTE'].copy()
                        p2 = routesites[-1]['PHEROIDS'].copy()
                        routes.append(p1)
                        pheromoneIDS.append(p2)
                        
                        #print(f"\n\n\n\nROUTES:{routes}\nPHEROIDS:{pheromoneIDS}")
                        #print("\n\n\n",routesites[-1]['PHEROIDS'].copy())
                        
                        route = []
                        feature_id_for_pheromone = []
                        route.append(backpoint)
                        site_tovisit.append(starting_id)
                        continue
                    elif backtracking and not site_tovisit:
                        #print("ROUTE COMPLETE")
                        route.remove(backpoint)
                        routes.append(route.copy())
                        pheromoneIDS.append(feature_id_for_pheromone.copy())
                        break
                        
                if current_node not in route:
                    route.append(current_node)
                
                #for traversal
                visibleRoads = []    #see the visible roads, this will be used to find roads which are intersected to the point   
                distancetoNextNode = []  
                routeToNextNode = []        #the route to the next node
                visibleNextNodes = []       #List of next nodes to visit
                pheroids = []               #Determines the Id's represented inside the pheromone layer
                pheroValues = []            #Values extracted from the ID's
                
                
                #this purpose was to iterate through the list of features and find the next node if the current node
                #intersected with the current feature line
                # Iterate through each road feature and check for intersection
                #these are feature lines which intersected to our point
                #print("DEBUG VISIBLE ROADS: ")
                visibleFeatureRoads = self.pointIntersect(current_node, None, self.roadLayer, flags=2)
                print(visibleFeatureRoads)
                #line of codes to be used that has more than 1 intersections
                if self.pointIntersect(current_node, None, self.roadLayer, flags=1) == 1 and not InitialDecision:
                    #print("CONTINUOUS")
                    ids = visibleFeatureRoads[0]
                    visibleFeature = self.roadLayer.getFeature(ids)
                    current_node = self.continousTraversal(
                        visibleFeature,
                        current_node,
                        route,
                        feature_id_for_pheromone,
                        target_site,
                        site_tovisit
                    )
                    continue
                else:
                    #print("intersect")
                    if InitialDecision:
                        InitialDecision = False
                    for ids in visibleFeatureRoads:
                        visibleFeature = self.roadLayer.getFeature(ids)
                        temproute = self.findNodeChoices(
                            visibleFeature,
                            current_node,
                            route,
                            target_site,
                            directionTol
                        )
                        #print("temproute NODECHOICES",temproute)
                        if temproute is None or not len(temproute):
                            continue
                        if type(temproute) is not list:
                            if temproute == current_node:
                                continue
                        if directionTol[0]:
                            for iroute in temproute:
                                if iroute is None:
                                    continue
                                visibleRoads.append(ids)
                                routeToNextNode.append(iroute)
                                visibleNextNodes.append(iroute[-1])
                                temppheroids = None
                                temppheroids = self.returnPheromoneEdgeID(iroute, flags = 1)
                                pheroids.append(temppheroids)
                                distancetoNextNode.append(self.returnHeuristicInfo(temppheroids,iroute[-1],site_tovisit,None))
                                for id in temppheroids:
                                    pheroValues.append(self.pheromoneLayer.getFeature(id).attribute("pheromone"))
                        elif type(temproute) == list:
                            if len(temproute)>1:
                                #print("more temproute")
                                for iroute in temproute:
                                    if iroute is None:
                                        continue
                                    elif iroute == current_node:
                                        continue
                                    visibleRoads.append(ids)
                                    visibleNextNodes.append(iroute)
                                    temppheroids = None
                                    temppheroids = self.returnPheromoneEdgeID([current_node,iroute], flags = 1)
                                    pheroids.append(temppheroids)
                                    distancetoNextNode.append(self.returnHeuristicInfo(temppheroids, iroute, site_tovisit,None))
                                    for id in temppheroids:
                                        pheroValues.append(self.pheromoneLayer.getFeature(id).attribute("pheromone"))
                                """print(f"VFR: {visibleRoads}",
                                    f"\nDNXTN: {distancetoNextNode}",
                                    f"\nVNN: {visibleNextNodes}",
                                    f"\nPHEROIDS: {pheroids}",
                                    f"\nPHEROVAL: {pheroValues}")"""
                            elif len(temproute):
                                #print("one temproute")
                                for iroute in temproute:
                                    if iroute is None:
                                        continue
                                    elif iroute == current_node:
                                        continue
                                    visibleRoads.append(ids)
                                    visibleNextNodes.append(iroute)
                                    temppheroids = None #TEMP
                                    temppheroids = self.returnPheromoneEdgeID([current_node,iroute], flags = 1)
                                    #print(f"ROUTE{iroute}\nIDS:{ids}, PHEROID:{temppheroids}")
                                    pheroids.append(temppheroids)
                                    distancetoNextNode.append(self.returnHeuristicInfo(temppheroids, iroute,site_tovisit,None))
                                    for id in temppheroids:
                                        pheroValues.append(self.pheromoneLayer.getFeature(id).attribute("pheromone"))
                                    """print(f"VFR: {visibleRoads}",
                                        f"\nDNXTN: {distancetoNextNode}",
                                        f"\nVNN: {visibleNextNodes}",
                                        f"\nPHEROIDS: {pheroids}",
                                        f"\nPHEROVAL: {pheroValues}")"""
                        else:
                            #print("solo bolo")
                            visibleRoads.append(ids)
                            visibleNextNodes.append(temproute)
                            temppheroids = None
                            temppheroids = self.returnPheromoneEdgeID([current_node,temproute], flags = 1)
                            pheroids.append(temppheroids)
                            distancetoNextNode.append(self.returnHeuristicInfo(temppheroids, temproute,site_tovisit,None))
                            for id in temppheroids:
                                pheroValues.append(self.pheromoneLayer.getFeature(id).attribute("pheromone"))
                print(f"VFR: {visibleRoads}\n RNN: {routeToNextNode}\n DNXTN: {distancetoNextNode} VNN: {visibleNextNodes}\nPHEROIDS: {pheroids} PHEROVAL: {pheroValues}")
                if visibleRoads:  # Check if the list is not empty
                    probabilities = np.array(pheroValues) ** self.alpha / np.array(distancetoNextNode) ** self.beta
                    probabilities /= probabilities.sum()
                    # Choose the next node based on probabilities
                    chosen_index = np.random.choice(len(visibleNextNodes), p=probabilities)
                    #print("CHOSEN INDEX: ",chosen_index)
                    chosen_node = visibleNextNodes[chosen_index]
                    current_node = chosen_node
                    feature_id_for_pheromone.extend(pheroids[chosen_index])
                    if directionTol[0]:
                        route.extend(routeToNextNode[chosen_index])
                        directionTol[0] = False
                    else:
                        route.append(current_node)
                    # Update current node and move to the next node
                else:
                    #print("pussy")
                    HARD_RESET += HARD_RESET * 1.15
                    #print(HARD_RESET)
                    if len(routesites)>0 and HARD_RESET > HARD_RESET_MAX and not backtracking:
                        print("\n\n\nHARD RESET with checkpoints")
                        temp = routesites.pop()
                        num_sites_visited.remove(temp['SITEID'])
                        site_tovisit.insert(temp["INDEX"],temp['SITEID'])
                        if temp['INDEX'] < iterate:
                            target_site = _CONST_SITES[temp['INDEX']]
                            iterate = temp['INDEX']
                        
                        if routesites:
                            temp = routesites[-1]   
                            route = temp["ROUTE"].copy()
                            feature_id_for_pheromone = temp["PHEROIDS"].copy()
                            current_node = route[-1]
                            
                        else:
                            print("big hard reset")
                            route, routesites = [], []
                            feature_id_for_pheromone = []
                            num_sites_visited=[]
                            site_tovisit = list(_CONST_SITES)
                            current_node = starting_node
                            iterate = 0
                            target_site = _CONST_SITES[0]
                        HARD_RESET = 1
                    elif len(routesites)>0 and not backtracking:
                        print("reset")
                        temp = routesites[-1]
                        route = temp["ROUTE"].copy()
                        feature_id_for_pheromone = temp["PHEROIDS"].copy()
                        current_node = route[-1]
                    else:
                        print("big hard reset")
                        route, routesites = [], []
                        feature_id_for_pheromone = []
                        num_sites_visited=[]
                        site_tovisit = list(_CONST_SITES)
                        current_node = starting_node
                        iterate = 0
                        target_site = _CONST_SITES[0]
                        HARD_RESET = 1
                self.iter_tour += 1
                #print(f"ROUTESITES:{routesites['SITEID']}\nPHEROIDS:{pheromoneIDS}")
                print(f"REMAIN SITES TO VISIT:{site_tovisit}TARGET SITE:{target_site}")
                
            # Store the route for the ant --
            route_line = QgsGeometry.fromPolylineXY(route)
            print(f"Ant {ant} is done traversing... ITERATION_TOUR: {self.iter_tour} ROUTE LENGTH {route_line.length()}")
            
            
            #THREAD RUN
            """result = ((route, distance), feature_id_for_pheromone)
            result_queue.put(result)"""
            
            
            return routes, pheromoneIDS
        def orderSiteNodes(self):
            siteNode_list = []
            for feature in self.siteNodes.getFeatures():
                if feature.attribute('start_node') == 'yes':
                    continue
                point = feature.geometry().asPoint()
                siteNode_list.append([
                    point,
                    feature.attribute("avg_volgen"),
                    feature.id()
                ])
            sorted_list = sorted(siteNode_list, key=lambda x: x[1], reverse=True)
            feature_list = []
            for fid in sorted_list:
                feature_list.append(fid[2])
            return feature_list
        def ACOAlgorithm(self, progress_callback=None):
            """The Main Part of the Algorithm, call to run ACO Algorithm
            
                Returns:
                    The best route
            
            """
            if self.pheromoneLayer is None:
                print("Initialize a pheromone Layer first before running ACO ALgorithm!")
                return None
            #VARIABLES FOR storing the best_route
            ITERATION_COUNTS = 0
            bestRoute = None
            bestDistance = float('inf')
            pheroIDs = []
            starting_id = None

            #Check start node
            site_nodes = self.orderSiteNodes()
            for node in self.siteNodes.getFeatures():
                if node['start_node'] == "yes":
                    starting_node = node.geometry().asPoint()
                    starting_id = node.id()
                    break

            #process of aco algorithm
            for iteration in range(self.num_iterations):
                ant_routes = []
                distance = float('inf')
                print(f"ITERATION {iteration+1} in process...")
                ITERATION_COUNTS = iteration+1
                
                ant_routes, feature_ids_pheromone=self.ant_traverse(starting_node, 1, site_nodes, starting_id, progress_callback)
                #print(ant_routes, feature_ids_pheromone)
                if self._stop_flag:
                    print('Interrupt')
                else:
                    print('All Ants have finished')
                # Find the best route among all ants
                """for (rd, fid_pher) in ant_routes:
                    min_distance = rd[1]
                    if min_distance < bestDistance:
                        bestDistance = min_distance
                        bestRoute = rd[0].copy()
                        pheroIDs = fid_pher.copy()"""

                # Update pheromone levels on road segments based on the feature id's
                #self.updatePheromones(pheroIDs,bestDistance)
            return ant_routes, 0.0
    
        def stop(self):
            self._stop_flag = True