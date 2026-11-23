# Usage Statistics

The charts below are populated live from the anonymised telemetry database. Data is only collected from users who have opted in. See [Telemetry](doc_telemetry.md) for details on what is collected and how to opt in or out. The page reads monthly totals only, never individual events. *Installs* is a lower bound: some older PRESTUS versions gave many installations the same ID.

<script src="https://cdn.jsdelivr.net/npm/@supabase/supabase-js@2"></script>
<script src="https://cdn.jsdelivr.net/npm/chart.js@4"></script>

<div id="telem-loading">Loading usage data…</div>
<div id="telem-error"></div>

<div class="telem-stats" id="telem-stats" style="display:none"></div>

<div class="telem-charts" id="telem-charts" style="display:none">
  <div class="telem-chart" style="grid-column:1/-1"><h4>Completed simulations per month</h4><canvas id="tc-timeline"></canvas></div>
  <div class="telem-chart"><h4>Execution platform</h4><canvas id="tc-platform"></canvas></div>
  <div class="telem-chart"><h4>Execution backend</h4><canvas id="tc-codetype"></canvas></div>
  <div class="telem-chart"><h4>Pipeline mode</h4><canvas id="tc-pipemode"></canvas></div>
  <div class="telem-chart"><h4>Simulation medium</h4><canvas id="tc-medium"></canvas></div>
  <div class="telem-chart"><h4>Mean run duration by medium</h4><canvas id="tc-duration"></canvas></div>
  <div class="telem-chart"><h4>Success rate by medium</h4><canvas id="tc-success-rate"></canvas></div>
</div>
