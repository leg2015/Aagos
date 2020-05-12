#ifndef AAGOS_WEB_D3_H
#define AAGOS_WEB_D3_H

/**
 * TODOs
 * - [ ] Allow world to be configured by web user
*/

#include "web/web.h"
#include "web/Input.h"

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
  UI::Document config_div;

  UI::Button run_toggle_but;
  UI::Button run_step_but;

  AagosPopulationVisualization pop_vis;

  void RedrawPopulation(bool update_data=true);
  void RedrawEnvironment();

public:
  AagosWebInterface(config_t & cfg)
    : AagosWorld(cfg),
      world_div("emp_world_view"),
      control_div("emp_controls_view"),
      config_div("emp_config_view"),
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
  run_toggle_but.SetAttr("class", "btn btn-block btn-lg btn-primary");

  // Step button setup.
  run_step_but = GetStepButton("run-step-button");
  run_step_but.SetAttr("class", "btn btn-block btn-lg btn-primary");

  // Add buttons to controls view.
  control_div << UI::Div("button-row").SetAttr("class", "row justify-content-md-center");

  control_div.Div("button-row")
    << UI::Div("step-col").SetAttr("class", "col-lg-auto p-2")
    << run_step_but;

  control_div.Div("button-row")
    << UI::Div("run-col").SetAttr("class", "col-lg-auto p-2")
    << run_toggle_but;

  control_div.Div("button-row")
    << UI::Div("render-frequency-col").SetAttr("class", "col-lg-auto p-2")
    << UI::Div("render-wrapper").SetAttr("class", "input-group input-group-lg")
    << UI::Div("render-input-prepend").SetAttr("class", "input-group-prepend")
    << UI::Div("render-input-prepend-text").SetAttr("class", "input-group-text")
    << "Render every";

  control_div.Div("render-wrapper")
    << UI::Input([](std::string in) { std::cout << "Change!" << std::endl; },
                 "number",
                 "",
                 "render-frequency")
        .SetAttr("class", "form-control")
        .SetCSS("min-width", "96px");

   control_div.Div("render-wrapper")
    << UI::Div("input-group-append").SetAttr("class", "input-group-append")
    << UI::Text().SetAttr("class", "input-group-text")
    << "th generation";

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

    // Configure
    pop_vis.Setup(*this);
    pop_vis.Start();

    emp::OnResize([this]() {
      std::cout << "Resize?" << std::endl;
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
    RedrawEnvironment();
  }
}

void AagosWebInterface::RedrawPopulation(bool update_data) {
  pop_vis.DrawPop(*this, update_data);
}

void AagosWebInterface::RedrawEnvironment() {
  pop_vis.DrawGradientEnv(*this, true);
}

#endif