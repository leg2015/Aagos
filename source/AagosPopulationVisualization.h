#ifndef AAGOS_WEB_POP_VIS_H
#define AAGOS_WEB_POP_VIS_H

#include "web/web.h"

#include "AagosWorld.h"
#include "AagosConfig.h"
#include "AagosOrg.h"

namespace UI = emp::web;

// NOTE - at the moment, this will explode if multiple instances of this object exist at once
class AagosPopulationVisualization  {
public:
  using world_t = AagosWorld;
  using org_t = typename AagosWorld::org_t;

  enum class POP_DRAW_MODE { FULL_POP=0, MAX_FIT };

protected:


  bool init=false;
  bool data_drawn=false;
  std::string element_id;

  UI::Document vis_div;

  size_t pop_org_height=30;     // height of organism in population view mode
  size_t max_fit_org_height=60; // Height of organism in max fit view mode
  size_t pop_view_max_height_px=500;

  POP_DRAW_MODE draw_mode=POP_DRAW_MODE::MAX_FIT;
  // POP_DRAW_MODE draw_mode=POP_DRAW_MODE::FULL_POP;

  void InitializeVariables(world_t & world) {
    // jswrap whatever it is that we need to jswrap

    // Initialize javascript object aagos_pop_vis
    EM_ASM({
      const elem_id = UTF8ToString($0);
      if (!emp.AagosPopVis) { emp.AagosPopVis = {}; }
      emp.AagosPopVis[elem_id] = {};
    }, element_id.c_str());
    // Initialize population data.
    InitializeGeneColorScale(world);
    InitializePopData();
    InitializeGradientEnvData();
  }

  void InitializeGeneColorScale(world_t & world) {
    EM_ASM({
      const elem_id = UTF8ToString($0);
      const num_genes = $1;
      // emp.AagosPopVis[elem_id]["gene_color_scale"] = d3.scaleSequential(d3.interpolateSinebow).domain([0, num_genes]);
      emp.AagosPopVis[elem_id]["gene_color_scale"] = d3.scaleSequential(d3.interpolateRainbow).domain([0, num_genes]);
      // emp.AagosPopVis[elem_id]["gene_color_scale"] = d3.scaleOrdinal().domain([0, num_genes]).range(d3["schemeAccent"]);
    }, element_id.c_str(),
       world.GetConfig().NUM_GENES());
  }

  void InitializePopData() {
    EM_ASM({
      var elem_id = UTF8ToString($0);
      emp.AagosPopVis[elem_id]["pop"] = [];
      emp.AagosPopVis[elem_id]["max_genome_size"] = 0;
      emp.AagosPopVis[elem_id]["most_fit_id"] = 0;
    }, element_id.c_str());
  }

  void InitializeGradientEnvData() {
    EM_ASM({
      var elem_id = UTF8ToString($0);
      emp.AagosPopVis[elem_id]["gradient_env"] = [];
      emp.AagosPopVis[elem_id]["num_gene_targets"] = 0;
      emp.AagosPopVis[elem_id]["gene_target_size"] = 0;
    }, element_id.c_str());
  }

  void UpdatePopDataFull(world_t & world) {
    if (!init || !world.IsSetup()) return;

    // info to update: bits, gene starts, num genes, gene size, num_bits
    EM_ASM({
      const elem_id = UTF8ToString($0);
      // Clear population information.
      emp.AagosPopVis[elem_id]['pop'] = [];
      emp.AagosPopVis[elem_id]['max_genome_size'] = 0;
      emp.AagosPopVis[elem_id]['most_fit_id'] = 0;

    }, element_id.c_str());

    size_t max_genome_size = 0;
    for (size_t org_id = 0; org_id < world.GetSize(); ++org_id) {
      org_t & org = world.GetOrg(org_id);
      const size_t num_bits = org.GetNumBits();
      const size_t gene_size = org.GetGeneSize();
      max_genome_size = num_bits > max_genome_size ? num_bits : max_genome_size;
      auto & gene_starts = org.GetGeneStarts();
      std::ostringstream stream;
      stream << org.GetBits();
      const std::string bits(stream.str());
      stream.str(std::string());
      for (size_t gene_id = 0; gene_id < gene_starts.size(); ++gene_id) {
        if (gene_id) stream << ",";
        stream << gene_starts[gene_id];
      }
      const std::string gene_starts_str(stream.str());
      EM_ASM({
        const elem_id = UTF8ToString($0);
        const org_id = $1;
        const bits = UTF8ToString($2).split('').map(Number);
        const gene_starts = UTF8ToString($3).split(',').map(Number);
        const gene_size = $4;
        const num_bits = bits.length;
        var gene_occupancy = new Map();
        var position_occupants = new Map();
        var gene_indicators = [];
        for (let gene_id = 0; gene_id < gene_starts.length; gene_id++) {
          gene_occupancy.set(gene_id, []);
          var start = gene_starts[gene_id];
          for (let k = 0; k < gene_size; k++) {
            const position = (start + k) % num_bits;
            gene_occupancy.get(gene_id).push(position);
            if (!(position_occupants.has(position))) {
              position_occupants.set(position, new Set());
            }
            gene_indicators.push({'gene_id': gene_id, 'pos': position, 'indicator_rank': position_occupants.get(position).size});
            position_occupants.get(position).add(gene_id);
          }
        }
        for (let i = 0; i < gene_indicators.length; i++) {
          gene_indicators[i]['rank']
        }
        emp.AagosPopVis[elem_id]['pop'].push({'org_id': org_id,
                                              'bits': bits,
                                              'gene_starts': gene_starts,
                                              'position_occupants': position_occupants,
                                              'gene_occupancy': gene_occupancy,
                                              'gene_indicators': gene_indicators});
      }, element_id.c_str(),       // 0
         org_id,                   // 1
         bits.c_str(),             // 2
         gene_starts_str.c_str(),  // 3
         gene_size);               // 4
    }

    EM_ASM({
      var elem_id = UTF8ToString($0);
      var most_fit_id = $1;
      var max_genome_size = $2;
      emp.AagosPopVis[elem_id]['most_fit_id'] = most_fit_id;
      emp.AagosPopVis[elem_id]['max_genome_size'] = max_genome_size;
    }, element_id.c_str(),   // 0
       world.GetMostFitID(), // 1
       max_genome_size);     // 2
  }

  // Do this separately instead of relying on D3 filters to save time!
  void UpdatePopDataMaxFit(world_t & world) {
    if (!init || !world.IsSetup()) return;

    const size_t most_fit_id = world.GetMostFitID();
    org_t & org = world.GetOrg(most_fit_id);
    const size_t genome_size = org.GetNumBits();
    const size_t gene_size = org.GetGeneSize();

    // Collect string representation of bits
    std::ostringstream stream;
    stream << org.GetBits();
    const std::string bits(stream.str());
    stream.str(std::string());
    // Collect string representation of gene starts
    auto & gene_starts = org.GetGeneStarts();
    for (size_t gene_id = 0; gene_id < gene_starts.size(); ++gene_id) {
      if (gene_id) stream << ",";
      stream << gene_starts[gene_id];
    }
    const std::string gene_starts_str(stream.str());

    // info to update: bits, gene starts, num genes, gene size, num_bits
    EM_ASM({
      const elem_id = UTF8ToString($0);
      const most_fit_id = $1;
      const genome_size = $2;
      const gene_starts = UTF8ToString($3).split(',').map(Number);
      const bits = UTF8ToString($4).split('').map(Number);
      const gene_size = $5;
      emp.AagosPopVis[elem_id]['most_fit_id'] = most_fit_id;
      emp.AagosPopVis[elem_id]['max_genome_size'] = genome_size;

      // Clear population information.
      emp.AagosPopVis[elem_id]['pop'] = [];

      var gene_occupancy = new Map();
      var position_occupants = new Map();
      var gene_indicators = [];
      for (let gene_id = 0; gene_id < gene_starts.length; gene_id++) {
        gene_occupancy.set(gene_id, []);
        var start = gene_starts[gene_id];
        for (let k = 0; k < gene_size; k++) {
          const position = (start + k) % genome_size;
          gene_occupancy.get(gene_id).push(position);
          if (!(position_occupants.has(position))) {
            position_occupants.set(position, new Set());
          }
          gene_indicators.push({'gene_id': gene_id, 'pos': position, 'indicator_rank': position_occupants.get(position).size});
          position_occupants.get(position).add(gene_id);
        }
      }
      for (let i = 0; i < gene_indicators.length; i++) {
        gene_indicators[i]['rank']
      }
      emp.AagosPopVis[elem_id]['pop'].push({'org_id': most_fit_id,
                                            'bits': bits,
                                            'gene_starts': gene_starts,
                                            'position_occupants': position_occupants,
                                            'gene_occupancy': gene_occupancy,
                                            'gene_indicators': gene_indicators});
    }, element_id.c_str(),      // 0
       most_fit_id,             // 1
       genome_size,             // 2
       gene_starts_str.c_str(), // 3
       bits.c_str(),            // 4
       gene_size);              // 5
  }

  void UpdateGradientEnvData(world_t & world) {
    // std::cout << "Update env data..." << std::endl;
    emp_assert(world.GetConfig().GRADIENT_MODEL());
    if (!init || !world.IsSetup()) return;
    const size_t num_gene_targets = world.GetConfig().NUM_GENES();
    const size_t gene_target_size = world.GetConfig().GENE_SIZE();
    // info to update: gradient env, num genes, gene sizes
    EM_ASM({
      const elem_id = UTF8ToString($0);
      const num_gene_targets = $1;
      const gene_target_size = $2;
      emp.AagosPopVis[elem_id]['gradient_env'] = [];
      emp.AagosPopVis[elem_id]['num_gene_targets'] = num_gene_targets;
      emp.AagosPopVis[elem_id]['gene_target_size'] = gene_target_size;
    }, element_id.c_str(), // 0
       num_gene_targets,   // 1
       gene_target_size);  // 2

    const auto & gene_targets = world.GetGradientFitnessModel().targets;
    for (size_t gene_id = 0; gene_id < num_gene_targets; ++gene_id) {
      std::ostringstream stream;
      stream << gene_targets[gene_id];
      const std::string gene_bits(stream.str());
      EM_ASM({
        const elem_id = UTF8ToString($0);
        const gene_bits = UTF8ToString($1).split('').map(Number);
        const gene_id = $2;
        emp.AagosPopVis[elem_id]["gradient_env"].push(({"bits":gene_bits,"id":gene_id}));
      }, element_id.c_str(),
         gene_bits.c_str(),
         gene_id);
    }
    // std::cout << "..done updating the environment" << std::endl;
  }

public:
  AagosPopulationVisualization(const std::string & id="pop-vis") : element_id(id), vis_div(id) { ; }

  void Setup(world_t & world) {
    std::cout << "Pop visualization setup" << std::endl;

    InitializeVariables(world);

    // Clean out contents of vis_div
    vis_div.Clear();

    // Add gene targets infrastructure
    vis_div << UI::Div(element_id + "-gene-targets-label-row").SetAttr("class", "row justify-content-center");
    vis_div.Div(element_id + "-gene-targets-label-row")
      << UI::Element("h5", element_id + "-gene-targets-label").SetAttr("class", "card-title")
      << "Gene Targets";

    vis_div << UI::Div(element_id + "-gene-targets-row").SetAttr("class", "row justify-content-center");
    vis_div.Div(element_id + "-gene-targets-row")
      << UI::Div(element_id + "-gene-targets-flex-row").SetAttr("class", "d-flex flex-row flex-wrap justify-content-center list-group-horizontal");

    for (size_t gene_id = 0; gene_id < world.GetConfig().NUM_GENES(); ++gene_id) {
      vis_div.Div(element_id + "-gene-targets-flex-row")
        << UI::Div(element_id + "-gene-target-" + emp::to_string(gene_id) + "-list-group-item").SetAttr("class", "list-group-item border border-dark rounded-0 m-1 p-0")
        << UI::Div(element_id + "-gene-target-" + emp::to_string(gene_id) + "-d-flex-container").SetAttr("class", "d-flex align-items-center h-100")
        << UI::Div(element_id + "-gene-target-" + emp::to_string(gene_id) + "-id-label").SetAttr("class", "d-flex align-items-center h-100 bg-dark text-light px-1")
        << gene_id;
      vis_div.Div(element_id + "-gene-target-" + emp::to_string(gene_id) + "-d-flex-container")
        << UI::Div(element_id + "-gene-target-" + emp::to_string(gene_id) + "-canvas-div");
      // NOTE - 'element_id + "-gene-target" + emp::to_string(gene_id) + "-canvas-div"' is where we'll
      //         put each gene target svg
    }

    vis_div << UI::Element("hr", element_id + "-gene-targets-hr");

    // Add population view infrastructure
    vis_div << UI::Div(element_id + "-population-label-row").SetAttr("class", "row justify-content-center");
    vis_div.Div(element_id + "-population-label-row")
      << UI::Element("h5", element_id + "-population-label").SetAttr("class", "card-title")
      << "Population";
    vis_div << UI::Div(element_id + "-population-canvas-row").SetAttr("class", "row");
    vis_div.Div(element_id + "-population-canvas-row")
      << UI::Div(element_id + "-population-canvas-div")
        .SetAttr("class", "AagosPopVis-canvas-div")
        .SetCSS("overflow-y", "scroll")
        .SetCSS("overflow-x", "scroll")
        .SetCSS("max-height", emp::to_string(pop_view_max_height_px) + "px");

    // Add SVG to #element_id
    EM_ASM({
      const elem_id = UTF8ToString($0);
      const num_genes = $1;

      // ---- Configure population view (div/svg/canvas/data-canvas) ----
      const pop_canvas_div_id = elem_id+"-population-canvas-div";
      emp.AagosPopVis[elem_id]["pop_canvas_div_id"] = pop_canvas_div_id;

      var pop_canvas_div = d3.select("#"+pop_canvas_div_id);
      pop_canvas_div.select("*").remove();

      var pop_svg = pop_canvas_div.append("svg")
                                  .attr("width",1)
                                  .attr("height",1)
                                  .attr("id","AagosPopVis-"+elem_id+"-pop-svg")
                                  .attr("class","AagosPopVis-svg");

      var pop_canvas = pop_svg.append("g")
                              .attr("id","AagosPopVis-"+elem_id+"-pop-canvas")
                              .attr("class","AagosPopVis-canvas");

      var pop_data_canvas = pop_canvas.append("g")
                                      .attr("id","AagosPopVis-"+elem_id+"-pop-data-canvas")
                                      .attr("class","AagosPopVis-canvas");

      // Setup gene target canvas/svg
      emp.AagosPopVis[elem_id]["gene_target_canvas_div_ids"] = [];
      for (let gene_id = 0; gene_id < num_genes; gene_id++) {
        const gene_target_canvas_div_id = elem_id + "-gene-target-" + gene_id + "-canvas-div";
        emp.AagosPopVis[elem_id]["gene_target_canvas_div_ids"].push(gene_target_canvas_div_id);

        var gene_target_canvas_div = d3.select("#"+gene_target_canvas_div_id);

        gene_target_canvas_div.select("*").remove(); // Clean out any old stuff.

        var gene_target_svg = gene_target_canvas_div.append("svg")
                                                    .attr("width", 1)
                                                    .attr("height", 1)
                                                    .attr("id", "AagosPopVis-"+elem_id+"-gene-target-"+gene_id+"-svg")
                                                    .attr("class", "AagosPopVis-svg");
        var gene_target_canvas = gene_target_svg.append("g")
                                                .attr("id", "AagosPopVis-"+elem_id+"-gene-target-"+gene_id+"-canvas");
        var gene_target_data_canvas = gene_target_canvas.append("g")
                                                        .attr("id","AagosPopVis-"+elem_id+"-gene-target-"+gene_id+"-data-canvas")
                                                        .attr("class","AagosPopVis-canvas");

      }

    }, element_id.c_str(), world.GetConfig().NUM_GENES());

    init = true;
  }

  void DrawGradientEnv(world_t & world, bool update_data=true) {
    if (update_data) UpdateGradientEnvData(world);
    EM_ASM({
      const elem_id = UTF8ToString($0);
      const bit_height = $1;

      var vis_info = emp.AagosPopVis[elem_id];
      var gene_colors = vis_info["gene_color_scale"];
      const num_gene_targets = vis_info["num_gene_targets"];
      const gene_target_size = vis_info["gene_target_size"];
      const max_genome_size = vis_info["max_genome_size"];

      const margins = ({top: 5, right: 5, bottom: 5, left: 5}); // todo - make class param

      if (!"org_bit_width" in vis_info) {
        vis_info["org_bit_width"] = 15; // todo - make class param
      }
      if (!"org_bit_height" in vis_info) {
        vis_info["org_bit_height"] = bit_height;
      }

      const org_bit_height = Math.min(bit_height, vis_info["org_bit_height"]);
      const org_bit_width = vis_info["org_bit_width"];

      const gene_target_width = gene_target_size * org_bit_width;
      const gene_target_height = org_bit_height;

      const canvas_width = gene_target_width + margins.left + margins.right;
      const canvas_height = gene_target_height + margins.top + margins.bottom;

      var gene_target_x_domain = ([0, gene_target_size]);
      var gene_target_x_range = ([0, canvas_width - margins.left - margins.right]);
      var x_scale = d3.scaleLinear().domain(gene_target_x_domain).range(gene_target_x_range);

      for (let gene_id = 0; gene_id < num_gene_targets; ++gene_id) {
        div_id = vis_info["gene_target_canvas_div_ids"][gene_id];
        svg_id = "AagosPopVis-"+elem_id+"-gene-target-"+gene_id+"-svg";
        canvas_id = "AagosPopVis-"+elem_id+"-gene-target-"+gene_id+"-canvas";
        data_canvas_id = "AagosPopVis-"+elem_id+"-gene-target-"+gene_id+"-data-canvas";

        // todo - update min width of div w/width of single gene + label!

        var svg = d3.select("#"+svg_id);
        svg.attr("width", canvas_width);
        svg.attr("height", canvas_height);

        var canvas = d3.select("#"+canvas_id);
        canvas.attr("transform", "translate(" + margins.left + "," + margins.top + ")");

        var data_canvas = d3.select("#"+data_canvas_id);
        data_canvas.selectAll("*").remove();

        var bits = data_canvas.selectAll("rect.bit").data(vis_info["gradient_env"][gene_id]["bits"]);
        bits.enter()
            .append("rect")
            .attr("class", "bit")
            .attr("id", function(bit, bit_i) {
              return "AagosPopVis-"+elem_id+"-gene-target-"+gene_id + "-bit-" + bit_i;
            })
            .attr("transform", function(bit, bit_i) {
              const x_trans = x_scale(bit_i);
              const y_trans = 0;
              return "translate(" + x_trans + "," + y_trans + ")";
            })
            .attr("height", org_bit_height)
            .attr("width", org_bit_width)
            .attr("fill", function(bit, bit_i) {
              return gene_colors(gene_id);
            })
            .attr("stroke", "gray");
        var bit_text = data_canvas.selectAll("text.bit").data(vis_info["gradient_env"][gene_id]["bits"]);
        bit_text.enter()
                .append("text")
                .attr("class", "bit-text")
                .attr("transform", function(bit, bit_i) {
                  const x_trans = x_scale(bit_i) + (org_bit_width / 2.0);
                  const y_trans = (org_bit_height / 2.0);
                  return "translate(" + x_trans + "," + y_trans + ")";
                })
                .attr("dominant-baseline", "middle")
                .attr("text-anchor", "middle")
                .text(function(bit, bit_i) {
                  return bit;
                });

      }

    }, element_id.c_str(),
       pop_org_height);
  }

  void DrawPop(world_t & world, bool update_data=true) {
    std::cout << "====== DRAWING POPULATION =========" << std::endl;
    if (update_data && draw_mode == POP_DRAW_MODE::FULL_POP) {
      UpdatePopDataFull(world);
    } else if (update_data && draw_mode == POP_DRAW_MODE::MAX_FIT) {
      UpdatePopDataMaxFit(world);
    }

    EM_ASM({
      // todo - move redundant stuff elsewhere/one-shot (no need to do multiple times)
      const elem_id = UTF8ToString($0);
      const org_height = $1;
      const pop_view_max_height_px = $2;

      var vis_info = emp.AagosPopVis[elem_id];
      const pop_size = vis_info["pop"].length;
      const max_genome_size = vis_info["max_genome_size"];
      const most_fit_id = vis_info["most_fit_id"];
      var gene_colors = vis_info["gene_color_scale"];

      const pop_div_id = emp.AagosPopVis[elem_id]["pop_canvas_div_id"];
      const pop_svg_id =  "#AagosPopVis-" + elem_id + "-pop-svg";
      const pop_canvas_id = "#AagosPopVis-" + elem_id + "-pop-canvas";
      const pop_data_canvas_id = "#AagosPopVis-" + elem_id + "-pop-data-canvas";

      const width = $('#' + elem_id).width(); // Width of surrounding div
      const margins = ({top:20, right:25, bottom:20, left:30}); // todo - make dynamic
      const height = (org_height * pop_size) + margins.top + margins.bottom;

      var canvas_height = height - margins.top - margins.bottom;

      console.log("height = " + height);

      min_bit_width = 15;
      min_canvas_width = (max_genome_size+1) * min_bit_width;

      var canvas_width = width - margins.left - margins.right;
      // todo - compute the min width based on a minimum bit size
      // If page is super small, allow left-right scrolling (don't shrink genomes too small)
      if (canvas_width < min_canvas_width) {
        canvas_width = min_canvas_width;
        // d3.select(pop_div_id).attr()
      }

      // Configure x/y range/domains (with a little bit of wiggle room)
      var pop_x_domain = ([0, max_genome_size]);
      var pop_x_range  = ([0, canvas_width]);
      var pop_y_domain = ([0, pop_size]);
      var pop_y_range  = ([0, canvas_height]);
      var pop_x_scale  = d3.scaleLinear().domain(pop_x_domain).range(pop_x_range);
      var pop_y_scale  = d3.scaleLinear().domain(pop_y_domain).range(pop_y_range);

      // Grab the pop svg to re-size appropriately!
      var svg = d3.select(pop_svg_id);
      svg.attr("width", canvas_width + margins.left + margins.right);
      svg.attr("height", height);

      const org_bit_width = pop_x_scale(1);
      const org_bit_height = pop_y_scale(0.9);
      vis_info["org_bit_width"] = org_bit_width;
      vis_info["org_bit_height"] = org_bit_height;

      // Grab the pop canvas inside of the svg, add transform to respect margins.
      var pop_canvas = d3.select(pop_canvas_id);
      pop_canvas.attr("transform", "translate(" + margins.left + "," + margins.top + ")");

      // --- Axes ---
      pop_canvas.selectAll(".axis").remove();
      if (pop_size > 1) {
        var y_axis = d3.axisLeft();
        y_axis.scale(pop_y_scale);
        pop_canvas.append("g").attr("class", "axis y_axis")
                              .attr("id", "AagosPopVis-"+elem_id+"-pop-y-axis")
                              .call(y_axis);
      }

      var x_axis = d3.axisTop();
      x_axis.scale(pop_x_scale);
      pop_canvas.append("g").attr("class", "axis x_axis")
                            .attr("id", "AagosPopVis-"+elem_id+"-pop-x-axis")
                            .call(x_axis);

      var axes = pop_canvas.selectAll(".axis");
      axes.selectAll("path").style({"fill": "none", "stroke": "black", "shape-rendering": "crispEdges"});
      axes.selectAll("text").style({"font-family": "sans-serif", "font-size": "10px"});


      var pop_data_canvas = d3.select(pop_data_canvas_id);
      pop_data_canvas.selectAll("*").remove();
      var organisms = pop_data_canvas.selectAll("g").data(vis_info["pop"]);
      organisms.enter()
               .append("g")
               .attr("class", "AagosPopVis-organism")
               .attr("id", function(org) {
                  return "AagosPopVis-"+elem_id+"-organism-"+org.org_id;
                })
               .attr("transform", function(org, org_i) {
                 const y_trans = pop_y_scale(org_i);
                 const x_trans = pop_x_scale(0);
                 return "translate(" + x_trans + "," + y_trans + ")";
                })
                .each(function(org) {
                      var bits = d3.select(this).selectAll("rect.bit").data(org["bits"]);
                      const rect_height = org_bit_height;
                      const rect_width = org_bit_width;
                      bits.enter()
                          .append("rect")
                          .attr("class", "bit")
                          .attr("transform", function(bit, bit_i) {
                            return "translate(" + pop_x_scale(bit_i) + ",0)";
                          })
                          .attr("height", rect_height)
                          .attr("width", rect_width)
                          .attr("fill", "white")
                          .attr("stroke", "gray");
                      var genes = d3.select(this).selectAll("rect.gene").data(org["gene_indicators"]);
                      genes.enter()
                           .append("rect")
                           .attr("class", "gene-indicator")
                           .attr("transform", function(indicator) {
                              const pos = indicator['pos'];
                              const num_occupants = org['position_occupants'].get(pos).size;
                              const rank = indicator['indicator_rank'];
                              const height = rect_height / num_occupants; // this should never be 0
                              const x_trans = pop_x_scale(pos);
                              const y_trans = height * rank;
                              return "translate(" + x_trans + "," + y_trans + ")";
                           })
                           .attr("height", function(indicator) {
                              const pos = indicator['pos'];
                              const num_occupants = org['position_occupants'].get(pos).size;
                              const height = rect_height / num_occupants;
                              return height;
                           })
                           .attr("width", rect_width)
                           .attr("fill", function(indicator) {
                             const gene_id = indicator['gene_id'];
                             return gene_colors(gene_id);
                           });
                      var bit_text = d3.select(this).selectAll("text").data(org["bits"]);
                      bit_text.enter()
                        .append("text")
                        .attr("class", "bit-text")
                        .attr("transform", function(bit, bit_i) {
                          const x_trans = pop_x_scale(bit_i) + (rect_width / 2.0);
                          const y_trans = (rect_height / 2.0);
                          return "translate(" + x_trans + "," + y_trans + ")";
                        })
                        .attr("dominant-baseline", "middle")
                        .attr("text-anchor", "middle")
                        .text(function(bit, bit_i) {
                          return bit;
                        });
                      });

    }, element_id.c_str(),
       draw_mode==POP_DRAW_MODE::FULL_POP ? pop_org_height : max_fit_org_height,
       pop_view_max_height_px);

    data_drawn=true;
  }
};

#endif