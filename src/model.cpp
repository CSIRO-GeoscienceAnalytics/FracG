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
	void BoundingBox_2D(Graph G, int nb_cells, int &p_tag, int &l_tag, float &lc)
	{
		std::vector<line_type> lines;

		for (auto Eg : make_iterator_range(edges(G)))
			lines.push_back(G[Eg].trace);  

		box AOI = ReturnAOI(lines);
	
		double min_x = bgm::get<bgm::min_corner, 0>(AOI) ;
		double min_y = bgm::get<bgm::min_corner, 1>(AOI) ;
		double max_x = bgm::get<bgm::max_corner, 0>(AOI) ;
		double max_y = bgm::get<bgm::max_corner, 1>(AOI) ;

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
	
std::vector<double> BoundingBox_3D(Graph G, int nb_cells, int &p_tag, int &l_tag, float &lc, double z)
	{


		
		
		double eps = 1e-2;
		std::vector<line_type> lines;

		for (auto Eg : make_iterator_range(edges(G)))
			lines.push_back(G[Eg].trace);  

		box AOI = ReturnAOI(lines);
	
		double min_x = bgm::get<bgm::min_corner, 0>(AOI) ;
		double min_y = bgm::get<bgm::min_corner, 1>(AOI) ;
		double max_x = bgm::get<bgm::max_corner, 0>(AOI) ;
		double max_y = bgm::get<bgm::max_corner, 1>(AOI) ;
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
		model::getEntitiesInBoundingBox(min_x-eps, min_y-eps, min_z-eps, max_x+eps, max_y+eps, max_z+eps, point_tags, 0);
		model::getEntitiesInBoundingBox(min_x-eps, min_y-eps, min_z-eps, max_x+eps, max_y+eps, max_z+eps, line_tags, 1);
		model::getEntitiesInBoundingBox(min_x-eps, min_y-eps, min_z-eps, max_x+eps, max_y+eps, max_z+eps, face_tags, 2);
		p_tag = point_tags.size();
		l_tag = line_tags.size();
		return(xyz);
	}
	
	
	
	void TagBoundaries_3D(std::vector<double> xyz)
	{
		float eps = 1e-2;
		gmsh::vectorpair volume_tags;
		double min_x = xyz.at(0);
		double max_x = xyz.at(1);
		double min_y = xyz.at(2);
		double max_y = xyz.at(3);
		double min_z = xyz.at(4);
		double max_z = xyz.at(5);
		
		model::getEntitiesInBoundingBox(min_x-eps, min_y-eps, min_z-eps, max_x+eps, max_y+eps, max_z+eps, volume_tags, 3);
		
		gmsh::vectorpair outDimTags;
		model::getBoundary(volume_tags, outDimTags, true, true, false);
		int nb_faces = outDimTags.size();

		std::vector<int> Left, Right, Front, Back, Top, Bottom;
		for (auto i : outDimTags) //loops through the surfaces (boundaries of volume)
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
					if (coord.at(c) == min_x)
						left +=1;
					if (coord.at(c) == max_x)
						right +=1;
					if (coord.at(c) == min_y)
						 front +=1;
					if (coord.at(c) == max_y)
						 back +=1;
					if (coord.at(c) == min_z)
						 bottom +=1;
					if (coord.at(c) == max_z)
						 top +=1;
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

		model::addPhysicalGroup(2, Left, (nb_faces + 1));
		model::setPhysicalName(2, (nb_faces + 1), "left");
		
		model::addPhysicalGroup(2, Right, (nb_faces + 2));
		model::setPhysicalName(2, (nb_faces + 2), "right");
		
		model::addPhysicalGroup(2, Front, (nb_faces + 3));
		model::setPhysicalName(2, (nb_faces + 3), "front");
		
		model::addPhysicalGroup(2, Back, (nb_faces + 4));
		model::setPhysicalName(2, (nb_faces + 4), "back");
		
		model::addPhysicalGroup(2, Top, (nb_faces + 5));
		model::setPhysicalName(2, (nb_faces + 5), "top");
		
		model::addPhysicalGroup(2, Bottom, (nb_faces + 6));
		model::setPhysicalName(2, (nb_faces + 6), "bottom");
		
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
	
	
	void MeshRefine_surface(std::vector<int> surf_tags, float lc, double gmsh_min_cl, double gmsh_min_dist, double gmsh_max_dist, std::vector<double> xyz )
	{
		double minCL, minDist, maxDist;


		//setting minimum characterisic length
		if (gmsh_min_cl == -1.0)
			minCL = lc/10; 
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
		
		std::vector<double> surf(surf_tags.begin(), surf_tags.end());
		
 
  		factory::synchronize();
  		

  		
		model::mesh::field::add("Distance", 1);
		//model::mesh::field::setNumber(1, "NNodesByEdge", 100);
		//model::mesh::field::setNumbers(1, "FacesList", surf);
		model::mesh::field::add("Threshold", 2);
		model::mesh::field::setNumber(2, "IField", 1);
		model::mesh::field::setNumber(2, "LcMin", minCL);
		model::mesh::field::setNumber(2, "LcMax", lc);
		model::mesh::field::setNumber(2, "DistMin", minDist);
		model::mesh::field::setNumber(2, "DistMax", maxDist);

		gmsh::model::mesh::field::add("Min", 3);
		gmsh::model::mesh::field::setNumbers(3, "FieldsList", {2});
		gmsh::model::mesh::field::setAsBackgroundMesh(3);
		
	
		gmsh::option::setNumber("Mesh.CharacteristicLengthExtendFromBoundary", 1);
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

	std::vector<int> EmbedLineaments_all(Graph G, std::vector< std::vector<std::pair<int, int>>> fused_lines, std::vector<int> intersec, int nb_bb_pts, std::string name)
	{
		std::ofstream ss_names, lowD_ss_name, interface_ss_name;
		std::string full_path  = FracG::AddPrefixSuffix(name, "", ".txt", true);
		std::string full_path2 = FracG::AddPrefixSuffix(name, "", "_lowD.txt", true);
		std::string full_path3 = FracG::AddPrefixSuffix(name, "", "_interface.txt", true);
		ss_names.open (full_path); 
		lowD_ss_name.open (full_path2); 
		interface_ss_name.open (full_path3); 
		factory::synchronize();

		std::vector<int> line_tags;
	   //get the tags of lineaments that are not the bounding box
	   for (int i = 0; i < num_edges(G); i++)
	   {
			std::vector<int>phys_group;
			std::vector<std::pair<int,int>> this_line_tags = fused_lines[i + nb_bb_pts];
			if (this_line_tags.empty())
				continue;

			for (auto ii : this_line_tags)
			{
				line_tags.push_back(ii.second);
				if (find(intersec.begin(), intersec.end(), ii.second) != intersec.end())
					phys_group.push_back(ii.second);
		}
		  model::addPhysicalGroup(1, phys_group, nb_bb_pts+i+1);
		  model::setPhysicalName(1, nb_bb_pts+i+1, "lineament_"+ to_string(i));
		  model::mesh::embed(1, phys_group, 2, 1);
		  ss_names << " lineament_"+ to_string(i);
		  lowD_ss_name << " lowerD_lineament_"+ to_string(i);
		  interface_ss_name << " interface_lineament_"+ to_string(i);
		}
		model::addPhysicalGroup(2, {1}, 1);
		model::setPhysicalName(2, 1, "host_rock");
		ss_names.close();
		lowD_ss_name.close(); 
		interface_ss_name.close();
		return(line_tags);
	}

	std::vector <box> CreateSamplingWindows(std::vector<line_type> faults, int nb_samples)
	{
		boost::random_device dev;
		boost::mt19937 rng(dev);

		std::vector <box> SamplingWindows;
		std::vector<line_type> lines;
		std::vector<double> length;

		BOOST_FOREACH(line_type l, faults)
		{
			lines.push_back(l);  
			length.push_back(bgm::length(l));
		}
		//find the shortest and longest lineament
		double min_len = floor(*min_element(length.begin(), length.end()) );
		double max_len = ceil(*max_element(length.begin(), length.end()) );

		float step = (max_len - min_len) / nb_samples;

		//create the bounding box for the entire set
		box AOI = ReturnAOI(lines);
		polygon_type t_AOI = ReturnTightAOI(lines);
		double min_x = bgm::get<bgm::min_corner, 0>(AOI);
		double min_y = bgm::get<bgm::min_corner, 1>(AOI);
		double max_x = bgm::get<bgm::max_corner, 0>(AOI);
		double max_y = bgm::get<bgm::max_corner, 1>(AOI);

		boost::random::uniform_real_distribution<double> rand_x(min_x,max_x);
		boost::random::uniform_real_distribution<double> rand_y(min_y,max_y);

		float size = min_len;

		for (int i = 0; i < nb_samples; i++)
		{
			bool found = false;
			do
			{
				point_type rand_p (rand_x(rng), rand_y(rng)); 
				point_type max_p ( (rand_p.x() + size), (rand_p.y() + size) );
				box sample_w{{rand_p.x(), rand_p.y()}, {max_p.x(), max_p.y()}}; 

				if (bgm::within(rand_p, t_AOI))
				{
					SamplingWindows.push_back(sample_w);
					found = true;
				}
			}while (!found);
		size += step;
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

	void AddLineament(line_type line, int source, int target, int &p_tag, int &l_tag, float lc)
	{
		int init_p_tag = p_tag + 1;
		std::vector<int> spline_points;

		bgm::unique(line);
		line_type simpl_line;
		bgm::simplify(line, simpl_line, 100);

		BOOST_FOREACH(point_type p, simpl_line)
			factory::addPoint(p.x(), p.y(), 0, ( lc ), ++p_tag);

		for( int i = init_p_tag; i < p_tag+1; i++ )
			spline_points.push_back( i );

		if (spline_points.size() < 2)
			std::cout << spline_points.size() << std::endl;
		factory::addSpline(spline_points, ++l_tag);
	}
	
	gmsh::vectorpair ExtrudeLines(gmsh::vectorpair line_tags, double x, double y, double z)
	{
		gmsh::vectorpair extruded_lines, new_entities;
		factory::extrude(line_tags, x,y,z, new_entities);
			
		for (auto newE : new_entities)
			if(newE.first == 2)
				extruded_lines.push_back(newE);
		
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
	
	std::vector<int> Fragment_and_cut(gmsh::vectorpair &toolDimTags, int bb_faces, float lc)
	{
		gmsh::vectorpair  dimTags;
		std::vector<std::pair<int, int> > ov, volume_tag, face_tag;
		std::vector<std::vector<std::pair<int, int> > > ovv;
		std::vector<int> v_tag;
		std::vector<std::vector<int> > f_tag;
		volume_tag.push_back(std::make_pair(3,1));
		
		factory::fragment(volume_tag, toolDimTags, ov, ovv, -1, true, true);
		factory::synchronize();
		GetTags(ovv, volume_tag, face_tag, v_tag, f_tag);
		
		model::addPhysicalGroup(3, v_tag, 1);
		model::setPhysicalName(3, 1, "host_rock");
		
		factory::intersect(volume_tag, face_tag, ov, ovv, -1, false, true);
		factory::synchronize();
		GetTags(ovv, volume_tag, face_tag, v_tag, f_tag);
		
	
		int i = 1;
		std::vector<int> surface_tag;
		for (auto Ftag : f_tag)
		{
			for (auto surface : Ftag)
				surface_tag.push_back(surface);
		}
		model::addPhysicalGroup(2, surface_tag, 10);
		model::setPhysicalName(2, 10, "surface_");
		return(surface_tag);
	}
	
	void WriteGmsh_2D(bool output, Graph G, int nb_cells, double gmsh_min_cl, double gmsh_min_dist, double gmsh_max_dist, std::string out_filename)
	{
		float lc;
		int nb_bb_pts;
		int p_tag = 0;
		int l_tag = 0;
		std::vector<int> intersec;
		std::vector< std::vector<std::pair<int, int>> > fused_lines;

		std::cout << "creating 2D mesh for lineament set" << std::endl;
		out_filename = FracG::AddPrefixSuffix(out_filename, "", ".msh");
		FracG::CreateDir(out_filename);
		gmsh::initialize();

		if (output)
			gmsh::option::setNumber("General.Terminal", 1);
		model::add("map2mesh");

		BoundingBox_2D(G, nb_cells, p_tag, l_tag, lc);

		nb_bb_pts = p_tag;
		for (auto Eg : make_iterator_range(edges(G)))
			AddLineament(G[Eg].trace, degree(source(Eg, G), G), degree(target(Eg, G),G), p_tag, l_tag, lc);
			
		MangeIntersections_bb(l_tag, nb_bb_pts, fused_lines);
		NameBoundingBox(nb_bb_pts, fused_lines, intersec);
		std::vector<int> line_tags = EmbedLineaments_all(G, fused_lines, intersec, nb_bb_pts, FracG::AddPrefixSuffix(out_filename, "", "_SideSet_names", true));
		MeshRefine_line(line_tags , lc, gmsh_min_cl, gmsh_min_dist, gmsh_max_dist);

		model::mesh::generate(2);
		gmsh::write(out_filename);
		if (output)
			gmsh::fltk::run();
		gmsh::finalize();
		std::cout << "Created msh-file " << out_filename << std::endl << std::endl;
	}
	
	
	void WriteGmsh_3D(bool output, Graph G, int nb_cells, double gmsh_min_cl, double gmsh_min_dist, double gmsh_max_dist, double z, std::string out_filename)
	{
		float lc;
		int nb_bb_pts;
		int p_tag = 0;
		int l_tag = 0;
		gmsh::vectorpair line_dim_tags, face_dim_tags;
		
		std::cout << "creating 3D mesh for lineament set" << std::endl;
		out_filename = FracG::AddPrefixSuffix(out_filename, "", ".msh");
		FracG::CreateDir(out_filename);
		gmsh::initialize();
		if (output)
			gmsh::option::setNumber("General.Terminal", 1);
		model::add("map2mesh");
		
		std::vector<double> xyz = BoundingBox_3D(G,  nb_cells, p_tag, l_tag, lc, z);
		
		nb_bb_pts = p_tag;
		for (auto Eg : make_iterator_range(edges(G)))
		{
			AddLineament(G[Eg].trace, degree(source(Eg, G), G), degree(target(Eg, G),G), p_tag, l_tag, lc);
			line_dim_tags.push_back(std::make_pair(1, l_tag));
		}
		factory::synchronize();
		
		face_dim_tags = ExtrudeLines(line_dim_tags, 0, 0, z);
		std::cout << "surfaces: " << face_dim_tags.size() << std::endl;
		
		factory::synchronize();
		std::vector<int> surf_tags = Fragment_and_cut(face_dim_tags, l_tag, lc);
		factory::synchronize();
	
		TagBoundaries_3D(xyz);
	
		MeshRefine_surface(surf_tags, lc, gmsh_min_cl, gmsh_min_dist, gmsh_max_dist, xyz);
		factory::synchronize();
		
		model::mesh::generate(3);
		gmsh::write(out_filename);
		
		std::cout << out_filename <<std::endl;
			gmsh::fltk::run();
		gmsh::finalize();
		std::cout << "Created msh-file " << out_filename << std::endl << std::endl;
	}
	
	
	
	void SampleNetwork_2D(bool output, VECTOR &lines, int nb_cells, int nb_samples, double map_distance_threshold, std::string filename)
	{
		std::vector<line_type> &faults = lines.data;
	// 	const char* dir = CreateDir(true);
		std::vector <box> sampling_windows;
		sampling_windows = CreateSamplingWindows(faults, nb_samples);

		std::cout << "Creating "<< nb_samples << " samples from lineament set" << std::endl;

		for (int w = 0; w < sampling_windows.size(); w++)
		{
			float lc;
			int nb_bb_pts;
			int p_tag = 0;
			int l_tag = 0;

			std::vector<int> intersec;
			std::vector< std::vector<std::pair<int, int>> > fused_lines;

			box AOI = sampling_windows.at(w);
			std::cout << AOI.min_corner().get<0>() << std::endl;
			Graph G = ReadVEC4MODEL(lines, AOI, map_distance_threshold);

			if (num_edges(G) > 0)
			{
				std::string output_filename = FracG::AddPrefixSuffix(filename, "", "_sample_" + to_string(w) +  ".msh");
				FracG::CreateDir(output_filename);
				gmsh::initialize();
				model::add(output_filename);
				if (output)
					gmsh::option::setNumber("General.Terminal", 1);

				double min_x = bgm::get<bgm::min_corner, 0>(AOI);
				double min_y = bgm::get<bgm::min_corner, 1>(AOI);
				double max_x = bgm::get<bgm::max_corner, 0>(AOI);
				double max_y = bgm::get<bgm::max_corner, 1>(AOI);

				lc = std::min((max_x-min_x)/nb_cells, (max_y-min_y)/nb_cells);

				factory::addPoint(min_x, min_y, 0, lc, ++p_tag);
				factory::addPoint(max_x, min_y, 0, lc, ++p_tag);
				factory::addPoint(max_x, max_y, 0, lc, ++p_tag);
				factory::addPoint(min_x, max_y, 0, lc, ++p_tag);

				for (int i = 1; i < p_tag; i++)
					factory::addLine(i, i+1, ++l_tag);
				factory::addLine(p_tag, 1, ++l_tag); 

				nb_bb_pts = p_tag;

				for (auto Eg : make_iterator_range(edges(G)))
					AddLineament(G[Eg].trace, degree(source(Eg, G), G), degree(target(Eg, G),G), p_tag, l_tag, lc);

				MangeIntersections_bb(l_tag,  nb_bb_pts, fused_lines);
				NameBoundingBox(nb_bb_pts, fused_lines, intersec);
				EmbedLineaments_all(G, fused_lines, intersec, nb_bb_pts, output_filename);

				factory::synchronize();
				model::mesh::generate(2);
				gmsh::write(output_filename);
				gmsh::finalize();
			}
		}
		std::cout << " done \n" << std::endl;
	}

}
