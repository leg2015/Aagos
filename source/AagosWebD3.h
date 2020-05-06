#ifndef AAGOS_WEB_D3_H
#define AAGOS_WEB_D3_H

/**
 * TODOs
 * - [ ] Allow world to be configured by web user

*/

#include "web/web.h"
#include "web/d3/dataset.h"

#include "AagosWorld.h"
#include "AagosConfig.h"
#include "AagosOrg.h"

namespace UI = emp::web;

double GetHTMLElementWidthByID(const std::string & id) {
  return EM_ASM_DOUBLE({
      var id = UTF8ToString($0);
      return $('#' + id).width();
    }, id.c_str());
}

double GetHTMLElementHeightByID(const std::string & id) {
  return EM_ASM_DOUBLE({
      var id = UTF8ToString($0);
      return $('#' + id).height();
    }, id.c_str());
}

// TODO - move to own file
// NOTE - at the moment, this will explode if multiple instances of this object exist at once
class AagosPopulationVisualization  {
public:
  using json_dataset_t = D3::JSONDataset;
  using world_t = AagosWorld;
  using org_t = typename AagosWorld::org_t;

protected:

  bool init=false;
  bool data_drawn=false;
  std::string element_id;

  void InitializeVariables() {
    // jswrap whatever it is that we need to jswrap

    // Initialize javascript object aagos_pop_vis
    EM_ASM({
      var elem_id = UTF8ToString($0);
      if (!emp.AagosPopVis) { emp.AagosPopVis = {}; }
    }, element_id.c_str());
    // Initialize population data.
    InitializePopData();
  }

  void InitializePopData() {
    std::cout << "Initialize pop data" << std::endl;
    EM_ASM({
      var elem_id = UTF8ToString($0);
      emp.AagosPopVis[elem_id] = {"pop":[]};
      emp.AagosPopVis[elem_id]["max_genome_size"] = 0;
      emp.AagosPopVis[elem_id]["most_fit_id"] = 0;
    }, element_id.c_str());
  }

  void UpdatePopData(world_t & world) {
    std::cout << "Update pop data..." << std::endl;
    if (!init) return;
    // info to update: bits, gene starts, num genes, gene size, num_bits
    EM_ASM({
      var elem_id = UTF8ToString($0);
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

public:
  AagosPopulationVisualization(const std::string & id="pop-vis") : element_id(id) { ; }

  void Setup() {
    std::cout << "Pop visualization setup" << std::endl;
    InitializeVariables();
    // Add SVG to #element_id
    EM_ASM({
      var elem_id = UTF8ToString($0);
      var svg = d3.select("#"+elem_id)
                  .append("svg").attr("width",1)
                                .attr("height",1)
                                .attr("id",elem_id+"-svg")
                                .attr("class","AagosPopVis-svg");
      svg.selectAll("*").remove();
      var canvas = svg.append("g").attr("id","AagosPopVis-"+elem_id+"-canvas")
                                  .attr("class","AagosPopVis-canvas");
      var env_canvas = canvas.append("g").attr("id","AagosPopVis-"+elem_id+"-env-canvas")
                                         .attr("class","AagosPopVis-canvas");
      var pop_canvas = canvas.append("g").attr("id","AagosPopVis-"+elem_id+"-pop-canvas")
                                         .attr("class","AagosPopVis-canvas");
      var pop_data_canvas = pop_canvas.append("g").attr("id","AagosPopVis-"+elem_id+"-pop-data-canvas")
                                                  .attr("class","AagosPopVis-canvas");
    }, element_id.c_str());
    // todo - setup auto-resizing
    init = true;
  }

  void Start() {
    // todo
  }

  void DrawPop(world_t & world, bool update_data=true) {
    std::cout << "draw pop" << std::endl;
    if (update_data) UpdatePopData(world); // todo - only do this if data is dirty
    EM_ASM({
      // todo - move redundant stuff elsewhere/one-shot (no need to do multiple times)
      var elem_id = UTF8ToString($0);
      var width = $('#' + elem_id).width();
      var height = 1000; // todo - make dynamic
      var margins = ({top:20,right:20,bottom:20,left:20}); // todo - make dynamic

      var vis = emp.AagosPopVis[elem_id];
      var max_genome_size = vis["max_genome_size"];
      var most_fit_id = vis["most_fit_id"];

      var vis_width = width;
      // todo - make all of this make sense for env + pop canvases!
      var canvas_width = vis_width - margins.left - margins.right;
      var canvas_height = height;

      var pop_x_domain = ([0, max_genome_size + 2]);
      var pop_x_range  = ([0, canvas_width]);
      var pop_y_domain = ([0, vis["pop"].length + 1]);
      var pop_y_range  = ([0, canvas_height]);
      var pop_x_scale  = d3.scaleLinear().domain(pop_x_domain).range(pop_x_range);
      var pop_y_scale  = d3.scaleLinear().domain(pop_y_domain).range(pop_y_range);

      var svg = d3.select("#"+elem_id+"-svg");
      svg.attr("width", canvas_width);
      svg.attr("height", canvas_height);

      var canvas = d3.select("#AagosPopVis-"+elem_id+"-canvas");
      canvas.attr("transform", "translate(" + margins.left + "," + margins.top + ")");
      var pop_canvas = d3.select("#AagosPopVis-"+elem_id+"-pop-canvas");

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

      var gene_colors = d3.scaleSequential(d3.interpolateSinebow).domain([0, 4]); // todo - get num genes

      var pop_data_canvas = d3.select("#AagosPopVis-"+elem_id+"-pop-data-canvas");
      pop_data_canvas.selectAll("*").remove();
      var organisms = pop_data_canvas.selectAll("g").data(vis["pop"]);
      organisms.enter()
               .append("g")
               .attr("class", "AagosPopVis-organism")
               .attr("id", function(d) {
                 return "AagosPopVis-"+elem_id+"-organism-"+d.org_id;
                })
               .attr("transform", function(d) {
                 var y_trans = pop_y_scale(d.org_id);
                 var x_trans = pop_x_scale(0);
                 return "translate(" + x_trans + "," + y_trans + ")";
                })
                .each(function(org, i) {
                      var bits = d3.select(this).selectAll("rect.bit").data(org["bits"]);
                      const rect_height = pop_y_scale(0.9);
                      const rect_width = pop_x_scale(1);
                      bits.enter()
                          .append("rect")
                          .attr("class", "bit")
                          .attr("transform", function(d, i) {
                            return "translate(" + pop_x_scale(i) + ",0)";
                          })
                          .attr("height", function(d, i) {
                            return rect_height;
                          })
                          .attr("width", function(d, i) {
                            return rect_width;
                          })
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
                           .attr("width", function(indicator) {
                             return rect_width;
                           })
                           .attr("fill", function(indicator) {
                             const gene_id = indicator['gene_id'];
                             return gene_colors(gene_id);
                           });
                      var bit_text = d3.select(this).selectAll("text").data(org["bits"]);
                      bit_text.enter()
                        .append("text")
                        .attr("class", "bit-text")
                        .attr("transform", function(d, i) {
                          var x_trans = pop_x_scale(i) + (rect_width/2.0);
                          var y_trans = (rect_height/2.0);
                          return "translate(" + x_trans + "," + y_trans + ")";
                        })
                        .attr("dominant-baseline", "middle")
                        .attr("text-anchor", "middle")
                        .text(function(d, i) {
                          return d;
                        });
                      });

    }, element_id.c_str()); // 0
    data_drawn=true;
  }
};

class AagosWebInterface : public UI::Animate, public AagosWorld {
public:
  // static constexpr double BIT_WIDTH
  static constexpr double BIT_HEIGHT = 20;
  static constexpr double MAX_GENE_TARGET_BIT_WIDTH = 10;
  static constexpr double GENE_TARGET_IDENTIFER_HEIGHT = 5;
  static constexpr double INDIV_VERT_MARGIN = 5;

  using config_t = AagosConfig;
  using canvas_draw_fun_t = std::function<void(void)>;

  using org_t = AagosOrg;
  using genome_t = AagosOrg::Genome;
  using phenotype_t = AagosOrg::Phenotype;

protected:

  UI::Document world_div;
  UI::Document control_div;

  UI::Button run_toggle_but;
  UI::Button run_step_but;

  emp::vector<UI::Canvas> env_canvases;
  AagosPopulationVisualization pop_vis;


  void RedrawPopulation(bool update_data=true);
  void RedrawEnvironment();

  void SetupEnvCanvasView();

  void ConfigEnvCanvasSize();
  void ConfigPopCanvasSize();

  void DrawPopCanvas_FullPop();
  void DrawPopCanvas_MaxFit();

  void DrawBit_SimpleBlackWhite();
  void DrawBit_SimpleText();

public:
  AagosWebInterface(config_t & cfg)
    : AagosWorld(cfg),
      world_div("emp_world_view"),
      control_div("emp_controls_view"),
      env_canvases(0),
      pop_vis("emp-pop-vis")
  {
    std::cout << "AagowWebInterface constructor!" << std::endl;
    SetupInterface();
  }

  /// Configure the interface
  void SetupInterface();

  /// Advance application by one frame.
  void DoFrame();
};

void AagosWebInterface::SetupInterface() {
  // ---- Setup interface control buttons ----
  // Run button setup.
  // Manually create run/pause button to have full control over button's callback.
  // NOTE - here's an example where Empirical's web tools are overly restrictive.
  run_toggle_but = UI::Button([this]() {
    ToggleActive();
    run_toggle_but.SetLabel(active ? "Stop" : "Run"); // TODO - use icons
    active ? run_step_but.SetDisabled(true) : run_step_but.SetDisabled(false);
  }, "Run", "run-toggle-button");
  run_toggle_but.SetAttr("class", "btn btn-primary m-1 ");

  // Step button setup.
  run_step_but = GetStepButton("run-step-button");
  run_step_but.SetAttr("class", "btn btn-primary m-1");

  // Add buttons to controls view.
  control_div
    << UI::Div("button-row").SetAttr("class", "row");

  control_div.Div("button-row")
    << UI::Div("step-col").SetAttr("class", "col")
    << run_step_but;

  control_div.Div("button-row")
    << UI::Div("run-col").SetAttr("class", "col")
    << run_toggle_but;

  // control_div.Div("button-row")
  //   << UI::Div("render-col").SetAttr("class", "col")
  //   <<

  // ---- Setup world view interface ----


  // Setup world view.
  // - environment -
  world_div << UI::Div("env-canvas-row").SetAttr("class", "row mb-2");
  world_div.Div("env-canvas-row")
    << UI::Div("env-canvas-col").SetAttr("class", "col");

  // - population -
  world_div << UI::Div("pop-canvas-row").SetAttr("class", "row");

  world_div.Div("pop-canvas-row")
    << UI::Div("pop-canvas-col").SetAttr("class", "col")
    << "<h4>Population</h4>"
    << UI::Div("emp-pop-vis");

  //   << pop_vis;
    // << pop_canvas;

  // pop_canvas_draw_fun = [this]() { DrawPopCanvas_FullPop(); };

  // Initial world configuration + pop canvas configuration.
  Setup(); // Call world setup


  // ---- Wire up event handlers ----

  UI::OnDocumentReady([this]() {
    std::cout << "-- OnDocumentReady (open) --"<<std::endl;
    SetupEnvCanvasView();
    ConfigEnvCanvasSize();
    ConfigPopCanvasSize();

    // Configure
    pop_vis.Setup();
    pop_vis.Start();

    emp::OnResize([this]() {
      std::cout << "Resize?" << std::endl;
      ConfigEnvCanvasSize();
      ConfigPopCanvasSize();
      RedrawPopulation(false);
      RedrawEnvironment();
    });

    RedrawPopulation(); // Finally, go ahead and draw initial population.
    RedrawEnvironment();
    std::cout << "-- OnDocumentReady (close) --"<<std::endl;
  });


}

void AagosWebInterface::DoFrame() {
  std::cout << "Frame!" << std::endl;
  RunStep();
  if (GetUpdate() % 50 == 0) {
    RedrawPopulation();
  }
}

void AagosWebInterface::RedrawPopulation(bool update_data) {
  // pop_canvas.Freeze();
  // pop_canvas_draw_fun();
  // pop_canvas.Activate();
  pop_vis.DrawPop(*this, update_data);
}

void AagosWebInterface::RedrawEnvironment() {
  emp_assert(config.GENE_SIZE() > 0);
  for (size_t gene_id = 0; gene_id < config.NUM_GENES(); ++gene_id) {
    auto & canvas = env_canvases[gene_id];
    auto & target = fitness_model_gradient->GetTarget(gene_id);
    canvas.Freeze();
    canvas.Clear("white");

    const double bit_width = canvas.GetWidth() / config.GENE_SIZE();
    const double bit_height = BIT_HEIGHT;
    const double gene_id_height = GENE_TARGET_IDENTIFER_HEIGHT;

    for (size_t bit_id = 0; bit_id < config.GENE_SIZE(); ++bit_id) {
      const double bit_x = bit_id * bit_width;
      const double bit_y = 0;
      const bool val = target.Get(bit_id);
      canvas.Rect(bit_x, bit_y, bit_width, bit_height, val ? "black" : "white", val ? "" : "black");
    }

    canvas.Activate();
  }
}

void AagosWebInterface::SetupEnvCanvasView() {
  world_div.Div("env-canvas-col").Clear();
  world_div.Div("env-canvas-col").SetAttr("class", "col");
  // Setup list group to hold gene target canvases.
  world_div.Div("env-canvas-col")
    << UI::Div("env-list-group").SetAttr("class", "list-group list-group-horizontal");
  world_div.Div("env-list-group")
    << UI::Div("gene-targets-label").SetAttr("class", "list-group-item bg-dark text-light")
    << "Gene Targets";
  // For each gene target, add canvas to env_canvases.
  emp_assert(config.GRADIENT_MODEL());
  env_canvases.clear();
  for (size_t i = 0; i < config.NUM_GENES(); ++i) {
    env_canvases.emplace_back(1, 1, "gene-target-canvas-" + emp::to_string(i));
    world_div.Div("env-list-group")
      << UI::Div("gene-target-" + emp::to_string(i)).SetAttr("class", "list-group-item")
      << env_canvases[i];
  }
}

void AagosWebInterface::ConfigEnvCanvasSize() {
  const double parent_width = GetHTMLElementWidthByID("emp_world_view");
  const double label_width = GetHTMLElementWidthByID("gene-targets-label"); // Reminder - this won't exist the first time this gets called during setup.
  const double workable_width = parent_width > label_width ? parent_width - label_width : parent_width;

  const double gene_canvas_width = emp::Min(workable_width / config.NUM_GENES(), MAX_GENE_TARGET_BIT_WIDTH * config.GENE_SIZE());
  const double gene_canvas_height = GENE_TARGET_IDENTIFER_HEIGHT + BIT_HEIGHT;

  for (size_t gene_id = 0; gene_id < env_canvases.size(); ++gene_id) {
    env_canvases[gene_id].SetSize(gene_canvas_width, gene_canvas_height);
  }
}

void AagosWebInterface::ConfigPopCanvasSize() {
  const size_t pop_size = GetSize();
  const double parent_width = GetHTMLElementWidthByID("emp_world_view");
  // pop_canvas.SetSize(parent_width, pop_size * (INDIV_VERT_MARGIN + BIT_HEIGHT));
}

#endif