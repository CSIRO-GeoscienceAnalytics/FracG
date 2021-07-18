#include "../include/model.h"
#include "../include/util.h"

namespace FracG
{
	namespace bgm = boost::geometry;
	
	const char* CreateDir(bool sampling)
	{
		const char* gmsh_dir = "./gmsh/";
		const char* gmsh_sam = "./gmsh/samples/";

		if (!sampling)
		{
			if(!opendir(gmsh_dir))
				mkdir(gmsh_dir, 0777);
			return (gmsh_dir);
		}
		else
		{
			if(!opendir(gmsh_sam))
				mkdir(gmsh_sam, 0777);
			return (gmsh_sam);
		}
	}
	
	
	/* create a bounding box of the area around the network.
	 * At the moment it is just a rectangular box derived from the envelope around all 
	 * lines */
	void BoundingBox_2D(Graph G, int nb_cells, double& x_m, double& y_m, int &p_tag, int &l_tag, float &lc, bool gmsh_in_meters)
	{
		std::vector<line_type> lines;

		double crop = 0;
		for (auto Eg : make_iterator_range(edges(G)))
			lines.push_back(G[Eg].trace);  

		box AOI = ReturnAOI(lines);
		if (gmsh_in_meters)
		{
			x_m = trunc(bgm::get<bgm::min_corner, 0>(AOI));
			y_m = trunc(bgm::get<bgm::min_corner, 1>(AOI));
		}
		else
		{
			x_m = 0;
			y_m = 0;
		}
		long double min_x  = bgm::get<bgm::min_corner, 0>(AOI) - x_m + crop ;
		long double min_y  = bgm::get<bgm::min_corner, 1>(AOI) - y_m + crop;
		long double max_x = bgm::get<bgm::max_corner, 0>(AOI)  - x_m - crop ;
		long double max_y = bgm::get<bgm::max_corner, 1>(AOI)  - y_m - crop;

		lc = std::min((max_x-min_x)/nb_cells, (max_y-min_y)/nb_cells);
		std::cout <<"init lc: " << lc<<std::endl;
		
		factory::addPoint(min_x, min_y, 0, lc, ++p_tag);
		factory::addPoint(max_x, min_y, 0, lc, ++p_tag);
		factory::addPoint(max_x, max_y, 0, lc, ++p_tag);
		factory::addPoint(min_x, max_y, 0, lc, ++p_tag);

		for (int i = 1; i < p_tag; i++)
			factory::addLine(i, i+1, ++l_tag);

		factory::addLine(p_tag, 1, ++l_tag); 
	}
	
	std::vector<double> BoundingBox_3D(Graph G, double& x_m, double& y_m, int nb_cells, int &p_tag, int &l_tag, float &lc, double z, bool gmsh_in_meters)
	{
		std::vector<line_type> lines;
		for (auto Eg : make_iterator_range(edges(G)))
			lines.push_back(G[Eg].trace);  
		double crop = 0;
		box AOI = ReturnAOI(lines);
	
		if (gmsh_in_meters)
		{
			x_m = trunc(bgm::get<bgm::min_corner, 0>(AOI));
			y_m = trunc(bgm::get<bgm::min_corner, 1>(AOI));
		}
		else
		{
			x_m = 0;
			y_m = 0;
		}
		
		long double min_x  = bgm::get<bgm::min_corner, 0>(AOI) - x_m + crop ;
		long double min_y  = bgm::get<bgm::min_corner, 1>(AOI) - y_m + crop;
		long double max_x = bgm::get<bgm::max_corner, 0>(AOI)  - x_m - crop ;
		long double max_y = bgm::get<bgm::max_corner, 1>(AOI)  - y_m - crop;

		double min_z = 0;
		double max_z = z;
		
		std::vector<double> xyz;
		xyz.push_back(min_x);
		xyz.push_back(max_x);
		xyz.push_back(min_y);
		xyz.push_back(max_y);
		xyz.push_back(min_z);
		xyz.push_back(max_z);
		
		lc = std::min((max_x-min_x)/nb_cells, (max_y-min_y)/nb_cells);
		std::cout <<"init lc: " << lc<<std::endl;
		
		int bb_tag = factory::addBox(min_x, min_y, min_z, (max_x - min_x), (max_y -min_y), z);
		factory::synchronize();
		
		gmsh::vectorpair point_tags, line_tags, face_tags;
		gmsh::model::getEntities(point_tags, 0);
		gmsh::model::getEntities(line_tags, 1);
		gmsh::model::getEntities(face_tags, 2);
		p_tag = point_tags.size();
		l_tag = line_tags.size();
		return(xyz);
	}

	void TagBoundaries_3D(std::vector<double> xyz, double gmsh_point_tol)
	{
		double dist = gmsh_point_tol;
		std::vector<int> Left, Right, Front, Back, Top, Bottom;
		gmsh::vectorpair volume_tags, surface_tags, outDimTags;;
		
		double min_x = xyz.at(0);
		double max_x = xyz.at(1);
		double min_y = xyz.at(2);
		double max_y = xyz.at(3);
		double min_z = xyz.at(4);
		double max_z = xyz.at(5);
		
		gmsh::model::getEntities(surface_tags, 2);
		gmsh::model::getEntities(volume_tags, 3);
		int nb_faces = surface_tags.size() ;

		model::getBoundary(volume_tags, outDimTags, true, true, false);
		
		for (auto i : outDimTags) //loops through the surfaces (boundaries of volumes)
		{
			gmsh::vectorpair in_tag, out_tag;
			in_tag.push_back(i);
			model::getBoundary(in_tag, out_tag, true, true, true); 
			
			int criteria;
			int left = 0, right = 0, front = 0, back = 0, top = 0, bottom = 0;
			for (auto j : out_tag) //loop through the points (boundaries of surfaces)
			{
				std::vector<double> parametricCoord;
				std::vector<double> coord;
				model::getValue(j.first, j.second, parametricCoord, coord);	
				criteria = out_tag.size();
				for (int c = 0; c < coord.size(); c++) //loop through the coodinate entries
				{
					if (abs(coord.at(c)-min_x) < dist)
						left +=1;
					if (abs(coord.at(c)-max_x) < dist)
						right +=1;
					if (abs(coord.at(c)-min_y) < dist)
						 front +=1;
					if (abs(coord.at(c)-max_y) < dist)
						 back +=1;
					if (abs(coord.at(c)-min_z) < dist)
						 top +=1;
					if (abs(coord.at(c)-max_z) < dist)
						 bottom +=1;
				}
			if(left == criteria)
				Left.push_back(i.second);
			if(right == criteria)
				Right.push_back(i.second);
			if(front == criteria)
				Front.push_back(i.second);
			if(back == criteria)
				Back.push_back(i.second);
			if(top == criteria)
				Top.push_back(i.second);
			if(bottom == criteria)
				Bottom.push_back(i.second);
			}
		}
	
			if(Left.size() != 0)
			{
				model::addPhysicalGroup(2, Left, (nb_faces + 100));
				model::setPhysicalName(2, (nb_faces + 100), "left");
			}
			if(Right.size() != 0)
			{
				model::addPhysicalGroup(2, Right, (nb_faces + 200));
				model::setPhysicalName(2, (nb_faces + 200), "right");
			}

			if(Front.size() != 0)
			{
				model::addPhysicalGroup(2, Front, (nb_faces + 300));
				model::setPhysicalName(2, (nb_faces + 300), "front");
			}
			
			if(Back.size() != 0)
			{
				model::addPhysicalGroup(2, Back, (nb_faces + 400));
				model::setPhysicalName(2, (nb_faces + 400), "back");
			}
			
			if(Top.size() != 0)
			{
				model::addPhysicalGroup(2, Top, (nb_faces + 500));
				model::setPhysicalName(2, (nb_faces + 500), "top");
			}
			
			if(Bottom.size() != 0)
			{
				model::addPhysicalGroup(2, Bottom, (nb_faces + 600));
				model::setPhysicalName(2, (nb_faces + 600), "bottom");
			}
		
	}

	void MeshRefine_line(std::vector<int> line_tags, float lc, double gmsh_min_cl, double gmsh_min_dist, double gmsh_max_dist)
	{
		double minCL, minDist, maxDist;
		std::vector<double> lines;
		std::vector<double> points;
		
		//setting minimum characterisic length
		if (gmsh_min_cl == -1.0)
			minCL = lc/10; 
		else
			minCL = gmsh_min_cl;

		//setting minimum distance for refinemnt around entity
		if (gmsh_min_dist == -1.0)
			minDist = lc/4; 
		else
			minDist = gmsh_min_dist;
		
		//setting maximum distance for refinemnt around entity
		if (gmsh_max_dist == -1.0)
			maxDist = lc/2;
		else
			maxDist = gmsh_max_dist;
		
		for (auto& i : line_tags)
			lines.push_back(i);		//need fixing

		factory::synchronize();
		model::mesh::field::add("Distance", 1);
		model::mesh::field::setNumber(1, "NNodesByEdge", 100);
		model::mesh::field::setNumbers(1, "EdgesList", lines);
		model::mesh::field::add("Threshold", 2);
		model::mesh::field::setNumber(2, "IField", 1);
		model::mesh::field::setNumber(2, "LcMin", minCL);
		model::mesh::field::setNumber(2, "LcMax", lc);
		model::mesh::field::setNumber(2, "DistMin", minDist);
		model::mesh::field::setNumber(2, "DistMax", maxDist);

		gmsh::model::mesh::field::add("Min", 3);
		gmsh::model::mesh::field::setNumbers(3, "FieldsList", {2});
		gmsh::model::mesh::field::setAsBackgroundMesh(3);

		gmsh::option::setNumber("Mesh.CharacteristicLengthExtendFromBoundary", 0);
		gmsh::option::setNumber("Mesh.CharacteristicLengthFromPoints", 0);
		gmsh::option::setNumber("Mesh.CharacteristicLengthFromCurvature", 0);
	}
	
	void MeshRefine_surface(std::vector<int> surf_tags, std::vector<double> edge_list, float lc, double gmsh_min_cl, double gmsh_min_dist, double gmsh_max_dist, std::vector<double> xyz )
	{
		double minCL, minDist, maxDist;

		//setting minimum characterisic length
		if (gmsh_min_cl == -1.0)
			minCL = lc/5; 
		else
			minCL = gmsh_min_cl;

		//setting minimum distance for refinemnt around entity
		if (gmsh_min_dist == -1.0)
			minDist = lc/8; 
		else
			minDist = gmsh_min_dist;
		
		//setting maximum distance for refinemnt around entity
		if (gmsh_max_dist == -1.0)
			maxDist = lc;
		else
			maxDist = gmsh_max_dist;
		
			std::vector<double> surf;
			for (auto i : surf_tags)
			{
				surf.push_back((double)i);
			}
			
  		factory::synchronize();
		model::mesh::field::add("Distance", 1);
		model::mesh::field::setNumber(1, "NNodesByEdge", 100);
		model::mesh::field::setNumbers(1, "EdgesList", edge_list);
		model::mesh::field::setNumbers(1, "FacesList", surf);
		model::mesh::field::add("Threshold", 2);
		model::mesh::field::setNumber(2, "IField", 1);
		model::mesh::field::setNumber(2, "LcMin", lc);
		model::mesh::field::setNumber(2, "LcMax", lc);
		model::mesh::field::setNumber(2, "DistMin", minDist);
		model::mesh::field::setNumber(2, "DistMax", maxDist);

		gmsh::model::mesh::field::add("Min", 3);
		gmsh::model::mesh::field::setNumbers(3, "FieldsList", {2});
		gmsh::model::mesh::field::setAsBackgroundMesh(3);

		gmsh::option::setNumber("Mesh.CharacteristicLengthExtendFromBoundary", 0);
		gmsh::option::setNumber("Mesh.CharacteristicLengthFromPoints", 0);
		gmsh::option::setNumber("Mesh.CharacteristicLengthFromCurvature", 0);
		
		gmsh::option::setNumber("Mesh.SmoothRatio", 5);
		gmsh::option::setNumber("Mesh.AnisoMax", 1000);
		gmsh::option::setNumber("Mesh.Algorithm", 6);

		gmsh::model::mesh::refine();
	}


	void NameBoundingBox(int nb_bb_pts, std::vector< std::vector<std::pair<int, int>> > fused_lines, std::vector<int>& intersec)
	{
		std::vector<int> fused_bb;
		std::vector<int> all_bb;

		std::vector<std::pair<int, int>> ov;
		std::vector <std::pair<int, int>> inside;
		std::vector<std::vector<std::pair<int, int>> > ovv;

		std::string boundary_tags[4] = {"bottom","right","top","left"};
		std::string phys_tag_name; 
		for (int i = 0; i < nb_bb_pts; i++)
		{
			phys_tag_name = "boundary_" + to_string(i);

			if (nb_bb_pts == 4)
				phys_tag_name = boundary_tags[i];

			for (auto ii : fused_lines[i])
				fused_bb.push_back(ii.second);
			all_bb.insert(all_bb.end(), fused_bb.begin(), fused_bb.end());

			model::addPhysicalGroup(1, fused_bb, i+1);
			model::setPhysicalName (1, {i+1}, phys_tag_name);
			fused_bb.clear();
		}

		factory::addCurveLoop(all_bb, 1);
		factory::addPlaneSurface({1}, 1);

		for (auto i : fused_lines)
		{
			for (auto ii : i)
				inside.push_back(std::make_pair(ii.first, ii.second));
		}
		factory::intersect(inside, {{2, 1}}, ov, ovv, -1, true, false);

		for (auto i : ov)
			intersec.push_back(i.second);
	}

	void MangeIntersections_bb(int &l_tag, int nb_bb_pts, std::vector< std::vector<std::pair<int, int>>>& fused_lines)
	{
	/*
	 * Mange intersections between lines and the bounding box
	 * and generate boundaries 
	*/
		std::vector<std::pair <int,int>> bb_lines;
		std::vector<std::pair <int,int>> l_lines;
		std::vector<std::pair <int,int>> ov;

		for( int i = 1; i < nb_bb_pts+1; i++ )
			bb_lines.push_back(std::make_pair(1, i));

		for( int i = nb_bb_pts+1; i < l_tag+1; i++ )
			l_lines.push_back(std::make_pair(1, i));

		factory::fuse(bb_lines, l_lines, ov, fused_lines, -1, true, true);
	}

	std::vector<int> EmbedLineaments_all(Graph G, std::vector< std::vector<std::pair<int, int>>> fused_lines, std::vector<int> intersec, int nb_bb_pts, std::string name, bool name_ss)
	{
			
		std::ofstream ss_names, lowD_ss_name, interface_ss_name, ss_properties;
		std::string full_path  = FracG::AddPrefixSuffix(name, "", ".txt", true);
		std::string full_path2 = FracG::AddPrefixSuffix(name, "", "_lowD.txt", true);
		std::string full_path3 = FracG::AddPrefixSuffix(name, "", "_interface.txt", true);
		std::string full_path4 = FracG::AddPrefixSuffix(name, "", "_propertiestxt", true);
		ss_names.open (full_path); 
		lowD_ss_name.open (full_path2); 
		interface_ss_name.open (full_path3); 
		ss_properties.open(full_path4);
		ss_properties <<"line_tag \t ID \t length \t parent_length \t angle \t set"<< std::endl;
		factory::synchronize();

		std::vector<int> line_tags;
		std::vector<int> all_lines_tags;
		
		std::vector<int> lines_tags_1;
		std::vector<int> lines_tags_2;
	   //get the tags of lineaments that are not the bounding box
		int i =0;
		for (auto Eg : make_iterator_range(edges(G)))
		{
			std::vector<int>phys_group;
			std::vector<std::pair<int,int>> this_line_tags = fused_lines[i + nb_bb_pts];
			if (this_line_tags.empty())
				continue;

			for (auto ii : this_line_tags)
			{
				line_tags.push_back(ii.second);
				phys_group.push_back(ii.second);
				all_lines_tags.push_back(ii.second);
			}
			
			if (name_ss)
			{
				model::addPhysicalGroup(1, phys_group, nb_bb_pts+i+1);
				model::mesh::embed(1, phys_group, 2, 1);
				ss_names << " line_"+ to_string(i);
				lowD_ss_name << " lowerD_line_"+ to_string(i);
				interface_ss_name << " interface_line_"+ to_string(i);
				ss_properties <<"line_"+ to_string(i) <<"\t"<< G[Eg].FaultNb <<"\t"<<  G[Eg].length <<"\t"<< G[Eg].fault_length << "\t"<< G[Eg].angle << "\t"<< G[Eg].set << std::endl;	
			}
		i++;
		}
		
		if(!name_ss)
		{
			model::addPhysicalGroup(1, all_lines_tags, nb_bb_pts+num_edges(G)+1);
			model::setPhysicalName(1, nb_bb_pts+num_edges(G)+1, "lines");
			model::mesh::embed(1, all_lines_tags, 2, 1);
			ss_names << " lines";
			lowD_ss_name << " lowerD_lines";
			interface_ss_name << " interface_lines";
			
		}
		
		model::addPhysicalGroup(2, {1}, 1);
		model::setPhysicalName(2, 1, "host_rock");
		ss_names.close();
		lowD_ss_name.close(); 
		interface_ss_name.close();
		return(line_tags);
	}

	std::vector <box> CreateSamplingWindows(Graph G, double gmsh_sw_size, int nb_samples)
	{
		boost::random_device dev;
		boost::mt19937 rng(dev);
		std::vector <box> SamplingWindows;
		std::vector<line_type> lines;
		std::vector<double> length;

		for (auto Eg : boost::make_iterator_range(edges(G))) 
		{
			lines.push_back(G[Eg].trace);
			length.push_back(bgm::length(G[Eg].trace));
		}

		//create the bounding box for the entire set
		box AOI = ReturnAOI(lines);
		polygon_type t_AOI = ReturnTightAOI(lines);
		double min_x = bgm::get<bgm::min_corner, 0>(AOI);
		double min_y = bgm::get<bgm::min_corner, 1>(AOI);
		double max_x = bgm::get<bgm::max_corner, 0>(AOI);
		double max_y = bgm::get<bgm::max_corner, 1>(AOI);
		
		point_type p1(min_x, min_y);
		point_type p2(max_x, max_y);
		double max_len = bgm::distance(p1,p2);
		double min_len = max_len / 4;
		
		boost::random::uniform_real_distribution<double> rand_x(min_x,max_x);
		boost::random::uniform_real_distribution<double> rand_y(min_y,max_y);
		boost::random::uniform_real_distribution<double> size(min_len, max_len);
		
		for (int i = 0; i < nb_samples; i++)
		{
			bool found = false;
			do
			{
				double d;
				if (gmsh_sw_size == -1)
					d = size(rng);
				else
					d = gmsh_sw_size;
				point_type rand_p (rand_x(rng), rand_y(rng)); 
				point_type max_p ( (rand_p.x() + d), (rand_p.y() + d) );
				box sample_w{{rand_p.x(), rand_p.y()}, {max_p.x(), max_p.y()}}; 

				if (!bgm::disjoint(sample_w, t_AOI))
				{
					SamplingWindows.push_back(sample_w);
					found = true;
				}
			}while (!found);
		}
		return(SamplingWindows);
	}

	p_tree BuildPointTree(Graph &G)
	{
		p_tree dist_tree;
		int index = 0;
		for (auto Eg : make_iterator_range(edges(G)))
		{
			BOOST_FOREACH(point_type p, G[Eg].trace)
				dist_tree.insert(std::make_pair(p, ++index));
		}
		return dist_tree;
	}
	
	void AddLineament(line_type line, double x_min, double y_min, int &p_tag, int &l_tag, float lc)
	{
		int init_p_tag = p_tag + 1;
		std::vector<int> spline_points;

		bgm::unique(line);
		line_type simpl_line;	
		bgm::simplify(line, simpl_line, 1);

		BOOST_FOREACH(point_type p, simpl_line)
		{
			factory::addPoint(p.x()-x_min, p.y()-y_min, 0, ( lc ), ++p_tag);
		}
		
		if (line.size() > 2)
		{
			for( int i = init_p_tag; i < p_tag+1; i++ )
				spline_points.push_back( i );

			if (spline_points.size() < 2)
				std::cout << spline_points.size() << std::endl;
			factory::addSpline(spline_points, ++l_tag);
		}
		if (line.size() == 2)
		{
			for( int i = init_p_tag; i < p_tag+1; i++ )
				spline_points.push_back( i );
			factory::addLine(spline_points[0],spline_points[1], ++l_tag);
		}
	}
	
	gmsh::vectorpair ExtrudeLines(FracG::AngleDistribution angle_dist, std::vector<std::pair<edge_type, int>> edge_map, Graph G, double z, std::vector<double> &edge_list)
	{
		gmsh::vectorpair extruded_lines, line_tags;
		std::cout << "lines: " << edge_map.size() << std::endl;
		for (auto edge_tag : edge_map)
		{
			gmsh::vectorpair cur_line, new_entities;
			cur_line.push_back(std::make_pair(1,edge_tag.second));
			//int set = CheckGaussians(angle_dist, G[edge_tag.first].angle);
			factory::extrude(cur_line, 0, 0, z, new_entities);
			
			for (auto newE : new_entities)
			{
				if(newE.first == 1)
					edge_list.push_back((double)newE.second);
				if(newE.first == 2)
					extruded_lines.push_back(newE);
			}
		}
		std::cout << "added surfaces: " << extruded_lines.size() << std::endl;
		return(extruded_lines);
	}
	
	void GetTags(std::vector<std::vector<std::pair<int, int> > > ovv, std::vector<std::pair<int, int> > &volume_tag, std::vector<std::pair<int, int> > &face_tag, std::vector<int> &v_tag, std::vector<std::vector<int> > &f_tag)
	{
		volume_tag.clear();
		face_tag.clear();
		v_tag.clear();
		f_tag.clear();
		for(auto ent : ovv)
		{
			std::vector<int> c_faces;
			for(auto i : ent)
			{
				if(i.first == 3)
				{
					volume_tag.push_back(std::make_pair(i.first, i.second));
					v_tag.push_back(i.second);
				}
				else if (i.first == 2)
				{
					face_tag.push_back(std::make_pair(i.first, i.second));
					c_faces.push_back(i.second);
				}
			}
			f_tag.push_back(c_faces);
		}
	}
	
	std::vector<int> Fragment_and_cut(gmsh::vectorpair &toolDimTags, int bb_faces, float lc, std::string name, bool name_ss)
	{
		std::ofstream ss_names, lowD_ss_name, interface_ss_name, ss_properties;
		std::string full_path  = FracG::AddPrefixSuffix(name, "", ".txt", true);
		std::string full_path2 = FracG::AddPrefixSuffix(name, "", "_lowD.txt", true);
		std::string full_path3 = FracG::AddPrefixSuffix(name, "", "_interface.txt", true);
		ss_names.open (full_path); 
		lowD_ss_name.open (full_path2); 
		interface_ss_name.open (full_path3); 
		
		gmsh::vectorpair  dimTags;
		std::vector<std::pair<int, int> > ov, volume_tag, face_tag;
		std::vector<std::vector<std::pair<int, int> > > ovv;
		std::vector<int> v_tag;
		std::vector<std::vector<int> > f_tag;

		gmsh::model::getEntities(volume_tag, 3);
		factory::fragment(volume_tag, toolDimTags, ov, ovv, -1, true, true);
		//factory::fragment(volume_tag, toolDimTags, ov, ovv);
		factory::synchronize();
		GetTags(ovv, volume_tag, face_tag, v_tag, f_tag);
		
		gmsh::model::getEntities(volume_tag, 3);
		int vn = 1;

		std::vector<int> volumes_in;
		for (auto v : volume_tag)
		{
			std::vector<int> cur_vol;
			cur_vol.push_back(v.second);
			volumes_in.push_back(v.second);
			model::addPhysicalGroup(3, cur_vol, vn);
			model::setPhysicalName(3, vn, "volume_"+ to_string(vn));
			vn++;
		}	
		
		//model::addPhysicalGroup(3, volumes_in, vn);
		//model::setPhysicalName(3, vn, "volume");
		
		factory::intersect(volume_tag, face_tag, ov, ovv, -1, false, true);
		factory::synchronize();
		GetTags(ovv, volume_tag, face_tag, v_tag, f_tag);
		
	
		int i = 1;
		std::vector<int> surface_tag;
		if (!name_ss)
		{
			for (auto Ftag : f_tag)
			{
				for (auto surface : Ftag)
					surface_tag.push_back(surface);
			}
			model::addPhysicalGroup(2, surface_tag, surface_tag.size());
			model::setPhysicalName(2, surface_tag.size(), "surfaces");
		}
		else
		{
			int tag = f_tag.size();
			for (auto Ftag : f_tag)
			{
				std::vector<int> surface_tag;
				for (auto surface : Ftag)
					surface_tag.push_back(surface);
					
				model::addPhysicalGroup(2, surface_tag, tag++);
				model::setPhysicalName(2, tag, "surface_"+ to_string(tag));
				
				ss_names << " surface_"+ to_string(tag);
				lowD_ss_name << " lowerD_surface_"+ to_string(tag);
				interface_ss_name << " interface_surface_"+ to_string(tag);
			}
		}
		return(surface_tag);
	}

	void WriteGmsh_2D(bool output, Graph G, int nb_cells, double gmsh_min_cl, double gmsh_min_dist, double gmsh_max_dist, bool gmsh_in_meters, bool name_ss, std::string out_filename)
	{
		float lc;
		int nb_bb_pts;
		int p_tag = 0;
		int l_tag = 0;
		double x_m,  y_m;
		std::vector<int> intersec;
		std::vector< std::vector<std::pair<int, int>> > fused_lines;

		std::cout << "creating 2D mesh for lineament set" << std::endl;
		out_filename = FracG::AddPrefixSuffix(out_filename, "", ".msh");
		FracG::CreateDir(out_filename);
		gmsh::initialize();

		if (output)
			gmsh::option::setNumber("General.Terminal", 1);
		model::add("map2mesh");

		gmsh::option::setNumber("General.NumThreads", 4);
		gmsh::option::setNumber("Geometry.Tolerance", 1e-15); 
		gmsh::option::setNumber("Geometry.MatchMeshTolerance", 1e-15); 
		BoundingBox_2D(G, nb_cells, x_m,  y_m, p_tag, l_tag, lc, gmsh_in_meters);

		nb_bb_pts = p_tag;
		for (auto Eg : make_iterator_range(edges(G)))
			AddLineament(G[Eg].trace, x_m,  y_m, p_tag, l_tag, lc);
			
		MangeIntersections_bb(l_tag, nb_bb_pts, fused_lines);
		
		NameBoundingBox(nb_bb_pts, fused_lines, intersec);
		std::vector<int> line_tags = EmbedLineaments_all(G, fused_lines, intersec, nb_bb_pts, FracG::AddPrefixSuffix(out_filename, "", "_SideSet_names", true), name_ss);
		MeshRefine_line(line_tags , lc, gmsh_min_cl, gmsh_min_dist, gmsh_max_dist);

		model::mesh::generate(2);
		gmsh::write(out_filename);
		if (output)
			gmsh::fltk::run();
		gmsh::finalize();
		std::cout << "Created msh-file " << out_filename << std::endl << std::endl;
	}
	
	void AddLayer_2D(gmsh::vectorpair& face_tags, std::vector<double> xyz, double lc, double depth, double tilt)
	{
		std::vector<std::pair<int, int> > volume_tags, surface_tags, new_face_tag;
		gmsh::model::getEntities(surface_tags, 2);
		gmsh::model::getEntities(volume_tags, 3);
		double x_min;
		if (depth == -700)
			x_min = xyz[0] + 4000;
		else 
			x_min = xyz[0] - 1000;
		double y_min = xyz[2] - 1000;
	
		double dx = (xyz[1] - xyz[0]) + 2000;
		double dy = (xyz[3] - xyz[2]) + 2000;
			
		double x_o = (x_min + (x_min+dx)) / 2;
		double y_o = (y_min + (y_min+dy)) / 2;
		std::cout << std::setprecision(10) << xyz[0] << " " << xyz[2] << " " << xyz[5] << std::endl;
		
		int surf = factory::addRectangle(x_min, y_min, depth, dx, dy, -1, 0.);
		face_tags.push_back(std::make_pair(2, surf));
		new_face_tag.push_back(std::make_pair(2, surf));
	
		std::cout << std::setprecision(10) << x_o << " " << y_o << " " << depth << std::endl;
			
		double dip = tilt * (3.14/180);
		factory::rotate(new_face_tag, x_o, y_o, depth, 1, 0, 0, dip);
		factory::synchronize();
		
		dip = 3 * (3.14/180);
		factory::rotate(new_face_tag, x_o, y_o, depth, 0, 1, 0, dip);
	}
	
	void WriteGmsh_3D(FracG::AngleDistribution angle_dist, bool output, Graph G, int nb_cells, double gmsh_min_cl, double gmsh_min_dist, double gmsh_max_dist, double z, double gmsh_point_tol, bool gmsh_in_meters, bool name_ss, std::string out_filename)
	{
		if (z > 0)
			z *= -1;
		float lc;
		int nb_bb_pts;
		int p_tag = 0;
		int l_tag = 0;
		double x_m,  y_m;
		
		std::vector<double>edge_list;
		gmsh::vectorpair line_dim_tags, face_dim_tags;
		std::vector<std::pair<edge_type, int>> edge_tag;
		
		std::cout << "creating 3D mesh for lineament set" << std::endl;
		out_filename = FracG::AddPrefixSuffix(out_filename, "", ".msh");
		FracG::CreateDir(out_filename);
		gmsh::initialize();
		if (output)
			gmsh::option::setNumber("General.Terminal", 1);
		model::add("map2mesh");
		
		gmsh::option::setNumber("General.NumThreads", 4);
		gmsh::option::setNumber("Geometry.Tolerance", 1e-15); 
		gmsh::option::setNumber("Geometry.MatchMeshTolerance", 1e-15); 
		
		std::vector<double> xyz = BoundingBox_3D(G, x_m, y_m, nb_cells, p_tag, l_tag, lc, z, gmsh_in_meters);

		nb_bb_pts = p_tag;
		for (auto Eg : make_iterator_range(edges(G)))
		{
			AddLineament(G[Eg].trace, x_m, y_m, p_tag, l_tag, lc);
			line_dim_tags.push_back(std::make_pair(1, l_tag));
			edge_tag.push_back(std::make_pair(Eg, l_tag));
		}
		factory::synchronize();
		
		face_dim_tags = ExtrudeLines(angle_dist, edge_tag, G, z, edge_list);
		std::cout << "surfaces: " << face_dim_tags.size() << std::endl;		
		factory::synchronize();
		
//Adding layeres--------------------------------------------------------
		 //AddLayer_2D(face_dim_tags, xyz, lc, -900, 2);
		//AddLayer_2D(face_dim_tags, xyz, lc, -1000, 5);
		//AddLayer_2D(face_dim_tags, xyz, lc, -700, 2);
		
		std::vector<int> surf_tags = Fragment_and_cut(face_dim_tags, l_tag, lc, out_filename, name_ss);
		factory::synchronize();
		TagBoundaries_3D(xyz, gmsh_point_tol);
		factory::synchronize();
			
		//MeshRefine_surface(surf_tags, edge_list, lc, gmsh_min_cl, gmsh_min_dist, gmsh_max_dist, xyz);
		factory::synchronize();
		
		model::mesh::generate(3);
		gmsh::write(out_filename);
		std::cout << out_filename <<std::endl;
		
		if (output)
			gmsh::fltk::run();
		gmsh::finalize();
		std::cout << "Created msh-file " << out_filename << std::endl << std::endl;
	}
				
	void SampleNetwork_2D(bool show_output, double gmsh_sw_size, Graph g, int nb_cells, int nb_samples, bool gmsh_in_meters, std::string filename)
	{
		if (nb_samples > 0)
		{
			std::vector <box> sampling_windows;
			sampling_windows = CreateSamplingWindows(g, gmsh_sw_size, nb_samples);
			std::cout << "Creating "<< nb_samples << " samples from lineament set" << std::endl;
			
			for (int w = 0; w < sampling_windows.size(); w++)
			{
				box AOI = sampling_windows.at(w);
				Graph G = Read4MODEL(g, AOI, 0.1);
				
				std::cout << "Min corner: " << std::setprecision(15) << AOI.min_corner().get<0>() <<" " << AOI.min_corner().get<1>() << std::endl;
				std::cout << "Max corner: " << std::setprecision(15) << AOI.max_corner().get<0>() <<" " << AOI.max_corner().get<1>() << std::endl;
				
				std::cout << " Contains "<< num_edges(G) << " edges" << std::endl;

				if (num_edges(G) > 0)
				{
					std::string output_filename = filename + "_sample_" + std::to_string(w);
					WriteGmsh_2D(show_output, G, 5, 50, 10, 20, gmsh_in_meters, true, output_filename );
				}
				else
					std::cout << "No features in sampling window  " << w << std::endl;
			}
			std::cout << " done \n" << std::endl;
		}
	}
}
