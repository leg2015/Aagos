#ifndef AAGOS_WEB_D3_H
#define AAGOS_WEB_D3_H

/**
 * TODOs
 * - [ ] Allow world to be configured by web user
*/

#include "web/web.h"

#include "AagosWorld.h"
#include "AagosConfig.h"
#include "AagosOrg.h"
#include "AagosPopulationVisualization.h"

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
    pop_vis.Setup(*this);
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
  pop_vis.DrawGradientEnv(*this, true);

  // emp_assert(config.GENE_SIZE() > 0);
  // for (size_t gene_id = 0; gene_id < config.NUM_GENES(); ++gene_id) {
  //   auto & canvas = env_canvases[gene_id];
  //   auto & target = fitness_model_gradient->GetTarget(gene_id);
  //   canvas.Freeze();
  //   canvas.Clear("white");

  //   const double bit_width = canvas.GetWidth() / config.GENE_SIZE();
  //   const double bit_height = BIT_HEIGHT;
  //   const double gene_id_height = GENE_TARGET_IDENTIFER_HEIGHT;

  //   for (size_t bit_id = 0; bit_id < config.GENE_SIZE(); ++bit_id) {
  //     const double bit_x = bit_id * bit_width;
  //     const double bit_y = 0;
  //     const bool val = target.Get(bit_id);
  //     canvas.Rect(bit_x, bit_y, bit_width, bit_height, val ? "black" : "white", val ? "" : "black");
  //   }

  //   canvas.Activate();
  // }
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