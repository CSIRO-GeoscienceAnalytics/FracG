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
	void BoundingBox_2d(Graph G, int nb_cells, int &p_tag, int &l_tag, float &lc)
	{
		std::vector<line_type> lines;

		for (auto Eg : make_iterator_range(edges(G)))
			lines.push_back(G[Eg].trace);  

		box AOI = ReturnAOI(lines);

		double min_x = bgm::get<bgm::min_corner, 0>(AOI);
		double min_y = bgm::get<bgm::min_corner, 1>(AOI);
		double max_x = bgm::get<bgm::max_corner, 0>(AOI);
		double max_y = bgm::get<bgm::max_corner, 1>(AOI);

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

	void MeshRefine_line(Graph G, float lc, int gmsh_min_cl, int gmsh_min_dist, int gmsh_max_dist)
	{
		std::vector<double> lines;
		std::vector<double> points;
		for (int i = 0; i < num_edges(G);i++)
			lines.push_back(double(i+4));		//need fixing

		int p_tag = 0;
		for (auto v : make_iterator_range(vertices(G)))
		{
			if(degree(v, G) > 1)
				points.push_back(double( p_tag+1 ));
			p_tag++;
		}

		gmsh::model::geo::synchronize();

		model::mesh::field::add("Distance", 1);
		// model::mesh::field::setNumber(1, "NNodesByEdge", 100);
		model::mesh::field::setNumbers(1, "EdgesList", lines);
		model::mesh::field::add("Threshold", 2);
		model::mesh::field::setNumber(2, "IField", 1);
		model::mesh::field::setNumber(2, "LcMin", lc/gmsh_min_cl);
		model::mesh::field::setNumber(2, "LcMax", lc);
		model::mesh::field::setNumber(2, "DistMin", lc/gmsh_min_dist);
		model::mesh::field::setNumber(2, "DistMax", lc/gmsh_max_dist);

		/*
		model::mesh::field::add("Distance", 10);
		model::mesh::field::setNumbers(10, "NodesList", points);
		model::mesh::field::add("Threshold", 20);
		model::mesh::field::setNumber(20, "IField", 10);
		model::mesh::field::setNumber(20, "LcMin", lc/10);
		model::mesh::field::setNumber(20, "LcMax", lc/2);
		model::mesh::field::setNumber(20, "DistMin", lc/4);
		model::mesh::field::setNumber(20, "DistMax", lc/2);
		*/

		gmsh::model::mesh::field::add("Min", 3);
		gmsh::model::mesh::field::setNumbers(3, "FieldsList", {2});
		gmsh::model::mesh::field::setAsBackgroundMesh(3);

		gmsh::option::setNumber("Mesh.CharacteristicLengthExtendFromBoundary", 0);
		gmsh::option::setNumber("Mesh.CharacteristicLengthFromPoints", 0);
		gmsh::option::setNumber("Mesh.CharacteristicLengthFromCurvature", 0);
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

	void EmbedLineaments_all(Graph G, std::vector< std::vector<std::pair<int, int>>> fused_lines, std::vector<int> intersec, int nb_bb_pts, std::string name)
	{
		std::ofstream ss_names, lowD_ss_name, interface_ss_name;
		std::string full_path  = FracG::AddPrefixSuffix(name, "", ".txt", true);
		std::string full_path2 = FracG::AddPrefixSuffix(name, "", "_lowD.txt", true);
		std::string full_path3 = FracG::AddPrefixSuffix(name, "", "_interface.txt", true);
		ss_names.open (full_path); 
		lowD_ss_name.open (full_path2); 
		interface_ss_name.open (full_path3); 
	   factory::synchronize();

	   //get the tags of lineaments that are not the bounding box
	   for (int i = 0; i < num_edges(G); i++)
	   {
			std::vector<int>phys_group;
			std::vector<std::pair<int,int>> this_line_tags = fused_lines[i + nb_bb_pts];
			if (this_line_tags.empty())
				continue;

			for (auto ii : this_line_tags)
				if (find(intersec.begin(), intersec.end(), ii.second) != intersec.end())
					phys_group.push_back(ii.second);

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
		int deg = 0; //TODO: deg needs to be set properly, this line only avoids undefined values
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

	void WriteGmsh_2D(bool output, Graph G, int nb_cells, int gmsh_min_cl, double gmsh_min_dist, double gmsh_max_dist, std::string out_filename)
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

		BoundingBox_2d(G, nb_cells, p_tag, l_tag, lc);

		nb_bb_pts = p_tag;
		for (auto Eg : make_iterator_range(edges(G)))
			AddLineament(G[Eg].trace, degree(source(Eg, G), G), degree(target(Eg, G),G), p_tag, l_tag, lc);

		MangeIntersections_bb(l_tag, nb_bb_pts, fused_lines);
		NameBoundingBox(nb_bb_pts, fused_lines, intersec);
		EmbedLineaments_all(G, fused_lines, intersec, nb_bb_pts, FracG::AddPrefixSuffix(out_filename, "", "_SideSet_names", true));
		MeshRefine_line(G, lc, gmsh_min_cl, gmsh_min_dist, gmsh_max_dist);

		model::mesh::generate(2);
		gmsh::write(out_filename);
		if (output)
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