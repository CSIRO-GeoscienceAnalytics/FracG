#include "../include/model.h"
#include "../include/util.h"

namespace FracG
{

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
		vector<line_type> lines;

		for (auto Eg : make_iterator_range(edges(G)))
			lines.push_back(G[Eg].trace);  

		box AOI = ReturnAOI(lines);

		double min_x = geometry::get<geometry::min_corner, 0>(AOI);
		double min_y = geometry::get<geometry::min_corner, 1>(AOI);
		double max_x = geometry::get<geometry::max_corner, 0>(AOI);
		double max_y = geometry::get<geometry::max_corner, 1>(AOI);

		lc = min((max_x-min_x)/nb_cells, (max_y-min_y)/nb_cells);
		cout <<"init lc: " << lc<<endl;

		factory::addPoint(min_x, min_y, 0, lc, ++p_tag);
		factory::addPoint(max_x, min_y, 0, lc, ++p_tag);
		factory::addPoint(max_x, max_y, 0, lc, ++p_tag);
		factory::addPoint(min_x, max_y, 0, lc, ++p_tag);

		for (int i = 1; i < p_tag; i++)
			factory::addLine(i, i+1, ++l_tag);

		factory::addLine(p_tag, 1, ++l_tag); 
	}

	void NameBoundingBox(int nb_bb_pts, vector< vector<pair<int, int>> > fused_lines, vector<int>& intersec)
	{
		vector<int> fused_bb;
		vector<int> all_bb;

		vector<pair<int, int>> ov;
		vector <pair<int, int>> inside;
		vector<vector<pair<int, int>> > ovv;

		string boundary_tags[4] = {"bottom","right","top","left"};
		string phys_tag_name; 
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
				inside.push_back(make_pair(ii.first, ii.second));
		}
		factory::intersect(inside, {{2, 1}}, ov, ovv, -1, true, false);

		for (auto i : ov)
			intersec.push_back(i.second);
	}

	void MangeIntersections_bb(int &l_tag, int nb_bb_pts, vector< vector<pair<int, int>>>& fused_lines)
	{
	/*
	 * Mange intersections between lines and the bounding box
	 * and generate boundaries 
	*/
		vector<pair <int,int>> bb_lines;
		vector<pair <int,int>> l_lines;
		vector<pair <int,int>> ov;

		for( int i = 1; i < nb_bb_pts+1; i++ )
			bb_lines.push_back(make_pair(1, i));

		for( int i = nb_bb_pts+1; i < l_tag+1; i++ )
			l_lines.push_back(make_pair(1, i));

		factory::fuse(bb_lines, l_lines, ov, fused_lines, -1, true, true);
	}

	void EmbedLineaments_all(Graph G, vector< vector<pair<int, int>>> fused_lines, vector<int> intersec, int nb_bb_pts, string name)
	{
		ofstream ss_names, lowD_ss_name, interface_ss_name;
		string full_path  = FracG::AddPrefixSuffix(name, "", ".txt", true);
		string full_path2 = FracG::AddPrefixSuffix(name, "", "_lowD.txt", true);
		string full_path3 = FracG::AddPrefixSuffix(name, "", "_interface.txt", true);
		ss_names.open (full_path); 
		lowD_ss_name.open (full_path2); 
		interface_ss_name.open (full_path3); 
	   factory::synchronize();

	   //get the tags of lineaments that are not the bounding box
	   for (int i = 0; i < num_edges(G); i++)
	   {
			vector<int>phys_group;
			vector<pair<int,int>> this_line_tags = fused_lines[i + nb_bb_pts];
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

	vector <box> CreateSamplingWindows(vector<line_type> faults, int nb_samples)
	{
		boost::random_device dev;
		boost::mt19937 rng(dev);

		vector <box> SamplingWindows;
		vector<line_type> lines;
		vector<double> length;

		BOOST_FOREACH(line_type l, faults)
		{
			lines.push_back(l);  
			length.push_back(geometry::length(l));
		}
		//find the shortest and longest lineament
		double min_len = floor(*min_element(length.begin(), length.end()) );
		double max_len = ceil(*max_element(length.begin(), length.end()) );

		float step = (max_len - min_len) / nb_samples;

		//create the bounding box for the entire set
		box AOI = ReturnAOI(lines);
		polygon_type t_AOI = ReturnTightAOI(lines);
		double min_x = geometry::get<geometry::min_corner, 0>(AOI);
		double min_y = geometry::get<geometry::min_corner, 1>(AOI);
		double max_x = geometry::get<geometry::max_corner, 0>(AOI);
		double max_y = geometry::get<geometry::max_corner, 1>(AOI);

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

				if (geometry::within(rand_p, t_AOI))
				{
					SamplingWindows.push_back(sample_w);
					found = true;
				}
			}while (!found);
		size += step;
		}
		return(SamplingWindows);
	}

	dist_tree BuildPointTree(Graph &G)
	{
		dist_tree dtree;
		int index = 0;
		for (auto Eg : make_iterator_range(edges(G)))
		{
			BOOST_FOREACH(point_type p, G[Eg].trace)
				dtree.insert(make_pair(p, ++index));
		}
		return dtree;
	}

	void AddLineament(line_type line, int source, int target, int &p_tag, int &l_tag, float lc, dist_tree &dtree)
	{
		int deg = 0; //TODO: deg needs to be set properly, this line only avoids undefined values
		int init_p_tag = p_tag + 1;
		vector<int> spline_points;

		geometry::unique(line);
		line_type simpl_line;

		boost::geometry::simplify(line, simpl_line, 100);

		BOOST_FOREACH(point_type p, simpl_line)
		{
	//find nearsest five neigbouring points and calculate average distance--
			vector<p_index> result;
			vector<double> distances;


	//mesh refinement-------------------------------------------------------
			dtree.query(geometry::index::nearest(p, 5), back_inserter(result));

			for (auto r = result.begin(); r < result.end(); r++)
				distances.push_back( geometry::distance(p, r[0].first));

			double average = accumulate( distances.begin(), distances.end(), 0.0)/distances.size();  
			double f = average/lc;
			//if ( f < 0.5)
			//	deg = average;
			//else
			//	deg = lc/floor(source + target) ; 

			if (geometry::equals(p, line.front()))
				deg /= source;

			if (geometry::equals(p, line.back()))
				deg /=  target;
	//----------------------------------------------------------------------
			factory::addPoint(p.x(), p.y(), 0, ( deg ), ++p_tag);
		}
		for( int i = init_p_tag; i < p_tag+1; i++ )
			spline_points.push_back( i );

		if (spline_points.size() < 2)
			cout << spline_points.size() << endl;
		factory::addSpline(spline_points, ++l_tag);
	}

	void WriteGmsh2D(dist_tree &dtree, bool output, Graph G, int nb_cells, string out_filename)
	{
	  cout << "creating 2D mesh for lineament set" << endl;
	  out_filename = FracG::AddPrefixSuffix(out_filename, "", ".msh");
	  FracG::CreateDir(out_filename);

	  float lc;
	  int nb_bb_pts;
	  int p_tag = 0;
	  int l_tag = 0;
	  vector<int> intersec;
	  vector< vector<pair<int, int>> > fused_lines;

	  gmsh::initialize();
	  if (output)
		gmsh::option::setNumber("General.Terminal", 1);
	  model::add("map2mesh");

	  BoundingBox_2d(G, nb_cells, p_tag, l_tag, lc);

	  nb_bb_pts = p_tag;

	  for (auto Eg : make_iterator_range(edges(G)))
		 AddLineament(G[Eg].trace, degree(source(Eg, G), G), degree(target(Eg, G),G), p_tag, l_tag, lc, dtree);

	  MangeIntersections_bb(l_tag, nb_bb_pts, fused_lines);
	  NameBoundingBox(nb_bb_pts, fused_lines, intersec);
	  EmbedLineaments_all(G, fused_lines, intersec, nb_bb_pts, FracG::AddPrefixSuffix(out_filename, "", "_SideSet_names", true));

	  factory::synchronize();
	  model::mesh::generate(2);
	  gmsh::write(out_filename);
	  if (output)
		gmsh::fltk::run();
	  gmsh::finalize();
	  cout << "Created msh-file " << out_filename << endl << endl;
	}

	void SampleNetwork2D(dist_tree &dtree, bool output, VECTOR &lines, int nb_cells, int nb_samples, double map_distance_threshold, string filename)
	{
		vector<line_type> &faults = lines.data;
	// 	const char* dir = CreateDir(true);
		vector <box> sampling_windows;
		sampling_windows = CreateSamplingWindows(faults, nb_samples);

		cout << "Creating "<< nb_samples << " samples from lineament set" << endl;

		for (int w = 0; w < sampling_windows.size(); w++)
		{
			float lc;
			int nb_bb_pts;
			int p_tag = 0;
			int l_tag = 0;

			vector<int> intersec;
			vector< vector<pair<int, int>> > fused_lines;

			box AOI = sampling_windows.at(w);
			cout << AOI.min_corner().get<0>() << endl;
			Graph G = ReadVEC4MODEL(lines, AOI, map_distance_threshold);

			if (num_edges(G) > 0)
			{
				string output_filename = FracG::AddPrefixSuffix(filename, "", "_sample_" + to_string(w) +  ".msh");
				FracG::CreateDir(output_filename);
				gmsh::initialize();
				model::add(output_filename);
				if (output)
					gmsh::option::setNumber("General.Terminal", 1);

				double min_x = geometry::get<geometry::min_corner, 0>(AOI);
				double min_y = geometry::get<geometry::min_corner, 1>(AOI);
				double max_x = geometry::get<geometry::max_corner, 0>(AOI);
				double max_y = geometry::get<geometry::max_corner, 1>(AOI);

				lc = min((max_x-min_x)/nb_cells, (max_y-min_y)/nb_cells);

				factory::addPoint(min_x, min_y, 0, lc, ++p_tag);
				factory::addPoint(max_x, min_y, 0, lc, ++p_tag);
				factory::addPoint(max_x, max_y, 0, lc, ++p_tag);
				factory::addPoint(min_x, max_y, 0, lc, ++p_tag);

				for (int i = 1; i < p_tag; i++)
					factory::addLine(i, i+1, ++l_tag);
				factory::addLine(p_tag, 1, ++l_tag); 

				nb_bb_pts = p_tag;

				for (auto Eg : make_iterator_range(edges(G)))
					AddLineament(G[Eg].trace, degree(source(Eg, G), G), degree(target(Eg, G),G), p_tag, l_tag, lc, dtree);

				MangeIntersections_bb(l_tag,  nb_bb_pts, fused_lines);
				NameBoundingBox(nb_bb_pts, fused_lines, intersec);
				EmbedLineaments_all(G, fused_lines, intersec, nb_bb_pts, output_filename);

				factory::synchronize();
				model::mesh::generate(2);
				gmsh::write(output_filename);
				gmsh::finalize();
			}
		}
		cout << " done \n" << endl;
	}

}