#ifndef AAGOS_WEB_POP_VIS_H
#define AAGOS_WEB_POP_VIS_H

#include "web/web.h"

#include "AagosWorld.h"
#include "AagosConfig.h"
#include "AagosOrg.h"

// namespace UI = emp::web;

// NOTE - at the moment, this will explode if multiple instances of this object exist at once
class AagosPopulationVisualization  {
public:
  using world_t = AagosWorld;
  using org_t = typename AagosWorld::org_t;

protected:

  bool init=false;
  bool data_drawn=false;
  std::string element_id;

  size_t org_height=30;
  size_t pop_view_max_height_px=500;

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
      emp.AagosPopVis[elem_id]["gene_color_scale"] = d3.scaleSequential(d3.interpolateSinebow).domain([0, num_genes]);
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

  void UpdatePopData(world_t & world) {
    std::cout << "Update pop data..." << std::endl;
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

  void UpdateGradientEnvData(world_t & world) {
    std::cout << "Update env data..." << std::endl;
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
    std::cout << "..done updating the environment" << std::endl;
  }

public:
  AagosPopulationVisualization(const std::string & id="pop-vis") : element_id(id) { ; }

  void Setup(world_t & world) {
    std::cout << "Pop visualization setup" << std::endl;
    InitializeVariables(world);
    // Add SVG to #element_id
    EM_ASM({
      var elem_id = UTF8ToString($0);
      var elem = d3.select("#"+elem_id);
      elem.selectAll("*").remove(); // Clean out root div.
      // ---- Configure environment view (div/svg/canvas/data-canvas) ----
      var env_div = elem.append("div")
                        .attr("id", "AagosPopVis-"+elem_id+"-env-canvas-div")
                        .style("overflow-x", "scroll");

      var env_svg = env_div.append("svg")
                           .attr("width",1)
                           .attr("height",1)
                           .attr("id","AagosPopVis-"+elem_id+"-env-svg")
                           .attr("class","AagosPopVis-svg");

      var env_canvas = env_svg.append("g")
                              .attr("id","AagosPopVis-"+elem_id+"-env-canvas")
                              .attr("class","AagosPopVis-canvas");

      var env_data_canvas = env_canvas.append("g")
                                      .attr("id","AagosPopVis-"+elem_id+"-env-data-canvas")
                                      .attr("class","AagosPopVis-canvas");

      // ---- Configure population view (div/svg/canvas/data-canvas) ----
      var pop_div = elem.append("div")
                        .attr("id", "AagosPopVis-"+elem_id+"-pop-canvas-div")
                        .style("overflow-y", "scroll");

      var pop_svg = pop_div.append("svg")
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
    }, element_id.c_str());
    // todo - setup auto-resizing
    init = true;
  }

  void Start() {
    // todo
  }

  void DrawGradientEnv(world_t & world, bool update_data=true) {
    std::cout << "draw env" << std::endl;
    if (update_data) UpdateGradientEnvData(world);
    EM_ASM({
      const elem_id = UTF8ToString($0);
      const org_height = $1;

      var vis_info = emp.AagosPopVis[elem_id];
      var gene_colors = vis_info["gene_color_scale"];
      const num_gene_targets = vis_info["num_gene_targets"];
      const gene_target_size = vis_info["gene_target_size"];
      const max_genome_size = vis_info["max_genome_size"];

      const env_div_id = "#AagosPopVis-" + elem_id + "-env-canvas-div";
      const env_svg_id =  "#AagosPopVis-" + elem_id + "-env-svg";
      const env_canvas_id = "#AagosPopVis-" + elem_id + "-env-canvas";
      const env_data_canvas_id = "#AagosPopVis-" + elem_id + "-env-data-canvas";

      var margins = ({top:20, right:20, bottom:20, left:20}); // todo - make dynamic
      const width = $('#' + elem_id).width(); // Width of surrounding div
      const canvas_width = width - margins.left - margins.right;

      // todo - handle case where canvas_width < gene_target_width
      const bit_width = width / (max_genome_size); // Approximately the same size as genome bits (a little bigger).
      const bit_height = org_height; // Approximately the same size as genome bits (a little bigger).
      const gene_target_h_padding = 10;
      const gene_target_v_padding = 5;

      const gene_target_width = (bit_width * (gene_target_size+2)); // spot for label + spacer
      var targets_per_row = Math.max(Math.floor(canvas_width / gene_target_width), 1);
      d3.select(env_div_id).style("width", width);
      const num_rows = Math.ceil(num_gene_targets / targets_per_row);

      const canvas_height = ((bit_height + gene_target_v_padding) * num_rows) + margins.top + margins.bottom;

      // Configure x/y range/domains (with a little bit of flex room?)
      var env_x_domain = ([0, targets_per_row]);
      var env_x_range  = ([0, canvas_width]);
      var env_y_domain = ([0, num_rows+1]);
      var env_y_range  = ([0, canvas_height]);
      var env_x_scale  = d3.scaleLinear().domain(env_x_domain).range(env_x_range);
      var env_y_scale  = d3.scaleLinear().domain(env_y_domain).range(env_y_range);

      var svg = d3.select(env_svg_id);
      svg.attr("width", canvas_width);
      svg.attr("height", canvas_height);

      var env_canvas = d3.select(env_canvas_id);
      env_canvas.attr("transform", "translate(" + margins.left + "," + margins.top + ")");

      env_canvas.selectAll(".env-canvas-title").remove();
      env_canvas.append("text")
                .attr("class", "env-canvas-title")
                .attr("transform", "translate(" + (canvas_width / 2) + ",-5)")
                .style("text-anchor", "middle")
                .text("Gene Targets");

      var env_data_canvas = d3.select(env_data_canvas_id);
      env_data_canvas.selectAll("*").remove();
      var gene_targets = env_data_canvas.selectAll("g").data(vis_info["gradient_env"]);
      gene_targets.enter()
                  .append("g")
                  .attr("class", "AagosPopVis-gene-target")
                  .attr("id", function(target) {
                    return "AagosPopVis-"+elem_id+"-gene-target-"+target.gene_id;
                  })
                  .attr("transform", function(target, gene_id) {
                    const row_pos = gene_id % targets_per_row;
                    const row_id = Math.floor(gene_id / targets_per_row);
                    const x_trans = env_x_scale(row_pos);
                    const y_trans = env_y_scale(row_id);
                    // console.log("targets_per_row=" + targets_per_row);
                    // console.log("row_pos=" + row_pos);
                    // console.log("row_id=" + row_id);
                    return "translate(" + x_trans + "," + y_trans + ")";
                  })
                  .each(function(target, gene_id) {
                    var bits = d3.select(this).selectAll("rect.bit").data(target["bits"]);
                    const rect_height = bit_height;
                    const rect_width = bit_width;
                    const row_pos = gene_id % targets_per_row;
                    const row_id = Math.floor(gene_id / targets_per_row);
                    bits.enter()
                        .append("rect")
                        .attr("class", "bit")
                        .attr("transform", function(bit, bit_i) {
                          const x_trans = (rect_width * bit_i); // +1.5 rw for label in front
                          const y_trans = 0;
                          return "translate(" + x_trans + "," + y_trans + ")";
                        })
                        .attr("height", rect_height)
                        .attr("width", rect_width)
                        .attr("fill", gene_colors(gene_id))
                        .attr('stroke', 'black');
                    var bit_text = d3.select(this).selectAll("text").data(target["bits"]);
                    bit_text.enter()
                            .append("text")
                            .attr("class", "bit-text")
                            .attr("transform", function(bit, bit_i) {
                              const x_trans = (rect_width * bit_i) + (rect_width / 2.0); // +1.5 rw for label in front
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
       org_height);
  }

  void DrawPop(world_t & world, bool update_data=true) {
    std::cout << "draw pop" << std::endl;
    if (update_data) UpdatePopData(world); // todo - only do this if data is dirty
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

      const pop_div_id = "#AagosPopVis-" + elem_id + "-pop-canvas-div";
      const pop_svg_id =  "#AagosPopVis-" + elem_id + "-pop-svg";
      const pop_canvas_id = "#AagosPopVis-" + elem_id + "-pop-canvas";
      const pop_data_canvas_id = "#AagosPopVis-" + elem_id + "-pop-data-canvas";

      const width = $('#' + elem_id).width(); // Width of surrounding div
      const height = org_height * pop_size;
      var margins = ({top:20, right:20, bottom:20, left:20}); // todo - make dynamic

      // If computed height is > minimum size, limit view height to maximum size.
      if (height > pop_view_max_height_px) {
        d3.select(pop_div_id).style("height", pop_view_max_height_px + "px");
      } else {
        d3.select(pop_div_id).style("height", "");
      }


      var canvas_width = width - margins.left - margins.right;
      var canvas_height = height - margins.top - margins.bottom;

      // If page is super small, allow left-right scrolling (don't shrink genomes too small)
      if (width < 350) {
        d3.select(pop_div_id).style("width", width);
        canvas_width = 350 - margins.left - margins.right;
      } else {
        d3.select(pop_div_id).style("width", "");
      }

      // Configure x/y range/domains (with a little bit of wiggle room)
      var pop_x_domain = ([0, max_genome_size + 2]);
      var pop_x_range  = ([0, canvas_width]);
      var pop_y_domain = ([0, pop_size + 1]);
      var pop_y_range  = ([0, canvas_height]);
      var pop_x_scale  = d3.scaleLinear().domain(pop_x_domain).range(pop_x_range);
      var pop_y_scale  = d3.scaleLinear().domain(pop_y_domain).range(pop_y_range);

      // Grab the pop svg to re-size appropriately!
      var svg = d3.select(pop_svg_id);
      svg.attr("width", canvas_width);
      svg.attr("height", canvas_height);

      // Grab the pop canvas inside of the svg, add transform to respect margins.
      var pop_canvas = d3.select(pop_canvas_id);
      pop_canvas.attr("transform", "translate(" + margins.left + "," + margins.top + ")");

      // Add label that will hangout in top margin
      pop_canvas.selectAll(".pop-canvas-title").remove();
      pop_canvas.append("text")
                .attr("class", "pop-canvas-title")
                .attr("transform", "translate(" + (canvas_width / 2) + ",-3)")
                .style("text-anchor", "middle")
                .text("Population");

      // --- Axes ---
      pop_canvas.selectAll(".y_axis").remove();
      var y_axis = d3.axisLeft();
      y_axis.scale(pop_y_scale);
      pop_canvas.append("g").attr("class", "axis y_axis")
                            .attr("id", "AagosPopVis-"+elem_id+"-pop-y-axis")
                            .call(y_axis);
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
               .attr("transform", function(org) {
                 const y_trans = pop_y_scale(org.org_id);
                 const x_trans = pop_x_scale(0);
                 return "translate(" + x_trans + "," + y_trans + ")";
                })
                .each(function(org) {
                      var bits = d3.select(this).selectAll("rect.bit").data(org["bits"]);
                      const rect_height = pop_y_scale(0.9);
                      const rect_width = pop_x_scale(1);
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
       org_height,
       pop_view_max_height_px);

    data_drawn=true;
  }
};

#endif