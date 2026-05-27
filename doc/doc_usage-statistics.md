# Usage Statistics

The charts below are populated live from the anonymised telemetry database. Data is only collected from users who have opted in. See [Telemetry](doc_telemetry.md) for details on what is collected and how to opt in or out.

<style>
.telem-stats { display:grid; grid-template-columns:repeat(auto-fit,minmax(140px,1fr)); gap:.75rem; margin:1.2rem 0; }
.telem-stat  { background:var(--md-code-bg-color,#f5f5f5); border-radius:8px; padding:1rem; text-align:center; }
.telem-stat .tv { font-size:1.8rem; font-weight:700; line-height:1.1; }
.telem-stat .tl { font-size:.72rem; color:var(--md-default-fg-color--light,#666); text-transform:uppercase; letter-spacing:.04em; margin-top:.2rem; }
.telem-charts { display:grid; grid-template-columns:repeat(auto-fit,minmax(340px,1fr)); gap:1rem; margin:1rem 0; }
.telem-chart  { background:var(--md-code-bg-color,#f5f5f5); border-radius:8px; padding:1rem; }
.telem-chart h4 { font-size:.75rem; font-weight:600; text-transform:uppercase; letter-spacing:.04em; color:var(--md-default-fg-color--light,#666); margin:0 0 .75rem; }
.telem-chart canvas { max-height:200px; }
#telem-error { display:none; color:#c0392b; font-size:.85rem; padding:.5rem 0; }
#telem-loading { color:var(--md-default-fg-color--light,#666); font-size:.85rem; }
</style>

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

<script>
(function () {
  const URL = 'https://lruelxwhkjibezkgipme.supabase.co';
  const KEY = 'sb_publishable_p18tafIOjIJuGwXb4szLSg_84vF3FL2';
  const PAL = ['#6366f1','#22d3ee','#f59e0b','#10b981','#f43f5e','#a855f7','#14b8a6','#fb923c','#84cc16','#e879f9'];

  // Auto-generated from git log — maps every 7-char commit hash to its tag family.
  // Regenerate with: git log --format="%h TAG" PREV_TAG..TAG for each tag pair.
  const HASH_INFO = {
    '00af2fd': 'v0.4.0', '00da4df': 'v0.6.0', '00e4861': 'v0.4.0', '013a6b6': 'v0.6.0', '01586cd': 'v0.6.0', '01869f6': 'v0.6.0',
    '018e78c': 'v0.1.0', '031d1fb': 'v0.3.0', '031f030': 'v0.3.0', '03409e5': 'v0.6.0', '034d4b3': 'v0.6.0', '03690fe': 'v0.6.0',
    '03afe3d': 'v0.6.0', '03dcacc': 'v0.4.0', '03f1b32': 'v0.4.0', '04183ac': 'v0.3.0', '043e9ee': 'v0.4.1', '045c6dc': 'v0.4.0',
    '04bfbf5': 'v0.6.0', '04c2bd6': 'v0.1.0', '04e62db': 'v0.5.0', '05146df': 'v0.2.1', '0528053': 'v0.6.0', '053771f': 'v0.6.0',
    '054181e': 'v0.1.0', '0546670': 'v0.4.0', '05f2f13': 'v0.2.0', '0662327': 'v0.3.0', '06847b2': 'v0.6.0', '06970d9': 'v0.1.0',
    '069c892': 'v0.6.1', '06ca159': 'v0.5.0', '06d1845': 'v0.4.0', '06e0a17': 'v0.4.2', '0704072': 'v0.6.0', '0732a6b': 'v0.5.0',
    '074629c': 'v0.4.2', '0772f85': 'v0.6.0', '07d93c8': 'v0.6.0', '07f3a23': 'v0.6.1', '0857dc2': 'v0.4.0', '0862dcf': 'v0.5.0',
    '08a20bb': 'v0.4.0', '08ee1d2': 'v0.2.0', '08fd90a': 'v0.6.0', '0928106': 'v0.6.0', '095d013': 'v0.2.0', '09f608c': 'v0.3.0',
    '0a3c25d': 'v0.4.2', '0a698cd': 'v0.6.0', '0a7a87f': 'v0.1.0', '0ab229c': 'v0.6.0', '0b42972': 'v0.6.0', '0b69017': 'v0.3.0',
    '0b89b97': 'v0.1.0', '0c1b5f6': 'v0.4.0', '0c2927b': 'v0.6.0', '0c62708': 'v0.6.0', '0c8cddc': 'v0.2.0', '0c973d1': 'v0.1.0',
    '0cc80ae': 'v0.5.0', '0cc8b0c': 'v0.5.0', '0d2ea98': 'v0.4.0', '0d54e8a': 'v0.4.0', '0d697f6': 'v0.1.0', '0d75ccb': 'v0.6.0',
    '0d76e58': 'v0.4.0', '0e9d222': 'v0.6.0', '0efb099': 'v0.3.0', '0f5a995': 'v0.5.0', '0f642ef': 'v0.6.0', '0f687c2': 'v0.6.0',
    '0f6a8aa': 'v0.1.0', '0f93dc2': 'v0.5.0', '0f97ee2': 'v0.6.1', '0fd714a': 'v0.6.0', '0fdb774': 'v0.6.0', '0ff68d0': 'v0.6.0',
    '10207fa': 'v0.5.0', '1082a43': 'v0.2.0', '109e2b4': 'v0.6.1', '10d34f1': 'v0.6.0', '10e82a1': 'v0.6.0', '110d0f3': 'v0.6.1',
    '1120090': 'v0.6.0', '1145913': 'v0.6.1', '1163742': 'v0.3.0', '116fe18': 'v0.4.0', '1171fdd': 'v0.6.0', '117336c': 'v0.1.0',
    '1174f80': 'v0.6.0', '11e3d0b': 'v0.6.0', '11fbeea': 'v0.1.0', '1208ac0': 'v0.6.0', '120aa39': 'v0.3.0', '1215130': 'v0.6.0',
    '12182b8': 'v0.3.0', '1224bff': 'v0.6.0', '126ae46': 'v0.3.0', '1286e50': 'v0.2.0', '12905a8': 'v0.6.0', '12d7165': 'v0.5.0',
    '13083f7': 'v0.1.0', '130e572': 'v0.6.0', '13a9069': 'v0.1.0', '13cec82': 'v0.6.0', '13e9461': 'v0.6.1', '1420b17': 'v0.2.0',
    '1448da5': 'v0.5.0', '144d49e': 'v0.6.1', '14676a5': 'v0.6.0', '14802b5': 'v0.5.0', '149e239': 'v0.6.0', '14b93ec': 'v0.1.0',
    '1558ec1': 'v0.3.0', '155a107': 'v0.4.0', '157f17e': 'v0.5.0', '159e4b8': 'v0.6.0', '15e6b56': 'v0.4.0', '1616f68': 'v0.4.0',
    '161a12a': 'v0.6.0', '1620676': 'v0.1.0', '1654eba': 'v0.6.0', '1682866': 'v0.6.0', '1692e14': 'v0.6.0', '16e04ca': 'v0.4.0',
    '17fde94': 'v0.6.0', '182aff2': 'v0.6.0', '184e37c': 'v0.6.1', '1931360': 'v0.6.0', '19683b7': 'v0.4.0', '19c1187': 'v0.6.0',
    '19db347': 'v0.3.0', '1a4ed69': 'v0.2.0', '1ad6c17': 'v0.4.0', '1b3aaf3': 'v0.6.0', '1b4afb4': 'v0.3.0', '1b5a2cc': 'v0.4.0',
    '1b95e30': 'v0.6.0', '1b9618c': 'v0.6.0', '1be0074': 'v0.6.0', '1c1443c': 'v0.4.0', '1c4410f': 'v0.5.0', '1c44185': 'v0.1.0',
    '1c558a1': 'v0.1.0', '1cb5aee': 'v0.6.0', '1d480bc': 'v0.5.0', '1dbac91': 'v0.6.1', '1dc62ee': 'v0.4.2', '1dcb42b': 'v0.5.0',
    '1dfae5c': 'v0.2.0', '1e3897c': 'v0.5.0', '1f0e3eb': 'v0.4.0', '1f70992': 'v0.3.0', '1f8b5c0': 'v0.6.0', '1fd43cd': 'v0.6.0',
    '1fe8995': 'v0.1.0', '1ff32bb': 'v0.4.2', '1ffd3dd': 'v0.4.2', '200bfc0': 'v0.6.0', '207567d': 'v0.6.0', '210df89': 'v0.6.0',
    '2115bca': 'v0.4.0', '2186ba0': 'v0.4.0', '225c44b': 'v0.1.0', '227d14e': 'v0.5.0', '22c2315': 'v0.1.0', '22d87b2': 'v0.6.0',
    '22dd02a': 'v0.6.0', '2348726': 'v0.5.0', '23b3e64': 'v0.6.0', '23d331e': 'v0.6.0', '242b54b': 'v0.4.0', '247d0f0': 'v0.6.1',
    '24c14b8': 'v0.6.0', '24c3cd1': 'v0.1.0', '24debe5': 'v0.4.0', '24e2461': 'v0.6.0', '25177ef': 'v0.6.0', '252712b': 'v0.6.1',
    '2555c83': 'v0.6.0', '2562651': 'v0.6.0', '25b9293': 'v0.6.0', '26241dc': 'v0.1.0', '2653900': 'v0.1.0', '26ab461': 'v0.5.0',
    '270e196': 'v0.6.1', '277e382': 'v0.6.1', '27c59a5': 'v0.6.0', '281e061': 'v0.5.0', '28daed1': 'v0.1.0', '28db2ff': 'v0.6.0',
    '291c5a6': 'v0.5.0', '294f797': 'v0.4.0', '29be0b4': 'v0.2.0', '29d89a0': 'v0.6.0', '2a93fc9': 'v0.6.0', '2b10558': 'v0.6.0',
    '2b4349d': 'v0.6.0', '2b99ce8': 'v0.6.0', '2c27f29': 'v0.3.0', '2c2b9c9': 'v0.6.0', '2c908aa': 'v0.6.0', '2d1ecae': 'v0.4.2',
    '2d7c014': 'v0.2.0', '2d7c21e': 'v0.2.0', '2db3126': 'v0.2.0', '2dba8ec': 'v0.6.0', '2e617dd': 'v0.6.0', '2e9bd07': 'v0.2.0',
    '2e9eadc': 'v0.6.0', '2ef1b7e': 'v0.6.0', '2f28df6': 'v0.4.2', '2f56459': 'v0.4.1', '2f84e92': 'v0.6.0', '2f97959': 'v0.1.0',
    '2ff8725': 'v0.6.0', '2ffdaf1': 'v0.6.0', '3025236': 'v0.6.0', '302f32a': 'v0.1.0', '30aec4e': 'v0.6.0', '30ed5f1': 'v0.5.0',
    '30f3c5c': 'v0.6.0', '312c62a': 'v0.3.0', '3139ce7': 'v0.4.0', '313eb4c': 'v0.5.0', '313fbff': 'v0.6.0', '3188796': 'v0.6.0',
    '3195d93': 'v0.3.0', '31ac292': 'v0.4.0', '31c0d64': 'v0.3.0', '31f166a': 'v0.4.0', '335a32a': 'v0.5.0', '34070dc': 'v0.4.0',
    '34315d2': 'v0.6.1', '345ecf6': 'v0.1.0', '3463e79': 'v0.4.0', '34ae7f0': 'v0.5.0', '34f9fc2': 'v0.6.0', '359ac7c': 'v0.6.0',
    '359ceff': 'v0.6.0', '35cebdb': 'v0.5.0', '35de8b1': 'v0.6.0', '361c869': 'v0.3.0', '365add2': 'v0.4.1', '365d3e3': 'v0.6.0',
    '367aff1': 'v0.4.0', '3686813': 'v0.6.1', '36b22e7': 'v0.5.0', '36b5980': 'v0.1.0', '36f734c': 'v0.6.0', '3710179': 'v0.4.0',
    '37210d8': 'v0.4.0', '3722539': 'v0.6.0', '3768eb4': 'v0.1.0', '378723e': 'v0.1.0', '386130e': 'v0.5.0', '38bfc46': 'v0.5.0',
    '38e0111': 'v0.6.0', '38f7088': 'v0.1.0', '38fe120': 'v0.4.1', '390397d': 'v0.5.0', '3969c4c': 'v0.6.1', '3978145': 'v0.4.0',
    '3989fe0': 'v0.2.0', '3993850': 'v0.1.0', '39c292f': 'v0.6.0', '39d61b0': 'v0.4.0', '39e453e': 'v0.3.0', '3a24937': 'v0.4.0',
    '3a82cca': 'v0.2.0', '3ae2c35': 'v0.4.2', '3af3d52': 'v0.5.0', '3b1a519': 'v0.5.0', '3b6821a': 'v0.6.0', '3b91ba7': 'v0.6.1',
    '3bb2dcc': 'v0.3.0', '3bd2316': 'v0.6.0', '3c15b88': 'v0.6.0', '3c89a18': 'v0.1.0', '3cac7f7': 'v0.6.0', '3ccec9f': 'v0.6.0',
    '3cef3b8': 'v0.5.0', '3d1df5d': 'v0.4.0', '3d59e6f': 'v0.2.0', '3dc6d65': 'v0.6.0', '3de52c4': 'v0.6.0', '3e429c1': 'v0.6.0',
    '3e4aa26': 'v0.6.0', '3e77571': 'v0.2.1', '3ee5256': 'v0.1.0', '3f099fa': 'v0.4.1', '3f3b24d': 'v0.6.0', '3fd5b93': 'v0.6.0',
    '40a01e6': 'v0.6.0', '40c5711': 'v0.6.0', '40d2024': 'v0.6.0', '40d251d': 'v0.6.0', '41820c4': 'v0.4.1', '41eac0c': 'v0.6.0',
    '41edd26': 'v0.6.0', '41f5542': 'v0.6.0', '4229129': 'v0.6.0', '4246df3': 'v0.4.0', '42aec82': 'v0.5.0', '42b51b1': 'v0.5.0',
    '436addb': 'v0.5.0', '4379632': 'v0.6.1', '43ea637': 'v0.6.0', '44090f2': 'v0.5.0', '441096e': 'v0.4.2', '4420015': 'v0.6.0',
    '4444e8e': 'v0.1.0', '44504f4': 'v0.5.0', '447e044': 'v0.6.0', '449bf12': 'v0.5.0', '44c71a8': 'v0.1.0', '450864b': 'v0.6.0',
    '457007d': 'v0.1.0', '45f1509': 'v0.6.0', '4682af3': 'v0.6.0', '471d449': 'v0.6.0', '47651f1': 'v0.6.0', '485c6b1': 'v0.3.0',
    '488bd14': 'v0.6.1', '48a496f': 'v0.5.0', '48b8f1c': 'v0.6.0', '48d4f89': 'v0.6.1', '491d602': 'v0.2.0', '49729d6': 'v0.3.0',
    '49be8b3': 'v0.4.0', '49feec7': 'v0.4.0', '4a44af3': 'v0.1.0', '4a5b79b': 'v0.1.0', '4a79a43': 'v0.6.0', '4aa2eb4': 'v0.6.0',
    '4acf1cf': 'v0.6.0', '4b3ebb2': 'v0.6.1', '4b59562': 'v0.5.0', '4b788d7': 'v0.6.1', '4b90619': 'v0.4.0', '4bbd13b': 'v0.4.0',
    '4c2cdfd': 'v0.4.0', '4c9cc44': 'v0.4.0', '4cb295a': 'v0.6.0', '4cf385c': 'v0.1.0', '4cfbbe0': 'v0.1.0', '4d1ad49': 'v0.6.0',
    '4d98718': 'v0.6.0', '4d9dcd4': 'v0.5.0', '4dec417': 'v0.6.0', '4e23050': 'v0.6.0', '4e30534': 'v0.6.1', '4e7ad86': 'v0.6.0',
    '4e939ac': 'v0.3.0', '4ea73f1': 'v0.6.0', '4f19af9': 'v0.6.1', '501ca12': 'v0.4.0', '502d9e8': 'v0.5.0', '505276b': 'v0.6.0',
    '5056506': 'v0.5.0', '509bea5': 'v0.6.1', '50dc480': 'v0.1.0', '512bbab': 'v0.6.0', '514692c': 'v0.1.0', '5148192': 'v0.5.0',
    '5193cf4': 'v0.6.0', '51a4935': 'v0.6.1', '52221eb': 'v0.6.0', '525f50c': 'v0.6.0', '52cb0b5': 'v0.6.0', '532c1d5': 'v0.2.0',
    '532ef3a': 'v0.4.0', '53a2ded': 'v0.6.0', '53f202b': 'v0.6.0', '540a1b9': 'v0.1.0', '548e2db': 'v0.3.0', '55036b1': 'v0.3.0',
    '5511162': 'v0.6.0', '5559d50': 'v0.6.0', '556daa2': 'v0.6.0', '55c2035': 'v0.3.0', '55c9cdb': 'v0.6.0', '561b5d6': 'v0.6.0',
    '569c8f0': 'v0.6.0', '56f6911': 'v0.4.0', '56ffa52': 'v0.3.0', '575f95c': 'v0.6.0', '5764dd6': 'v0.5.0', '57b2be1': 'v0.5.0',
    '57c8798': 'v0.6.0', '57d1a16': 'v0.6.0', '58241e2': 'v0.4.0', '5849299': 'v0.5.0', '58a1a8b': 'v0.5.0', '58a6388': 'v0.4.0',
    '58c5a60': 'v0.6.1', '59155af': 'v0.6.0', '5949d3b': 'v0.6.0', '59c8c1a': 'v0.3.0', '59e0f16': 'v0.5.0', '59f7143': 'v0.5.0',
    '5a1e095': 'v0.3.0', '5a4506f': 'v0.6.0', '5a5c48e': 'v0.6.0', '5a88d81': 'v0.6.0', '5a8a2e9': 'v0.2.0', '5ac113f': 'v0.4.0',
    '5b0a80b': 'v0.4.0', '5b43dfb': 'v0.5.0', '5b4596e': 'v0.2.0', '5b45e88': 'v0.6.0', '5b83233': 'v0.4.2', '5b83e13': 'v0.6.0',
    '5bc6584': 'v0.4.0', '5bdfcff': 'v0.6.0', '5bff713': 'v0.5.0', '5c05cc7': 'v0.6.0', '5c4a8b9': 'v0.4.0', '5cc9d4f': 'v0.6.0',
    '5ce42a7': 'v0.5.0', '5cfbe37': 'v0.2.0', '5d22507': 'v0.3.0', '5d3110c': 'v0.6.1', '5d512b8': 'v0.6.0', '5e1dd92': 'v0.6.1',
    '5eb5730': 'v0.6.0', '5eba932': 'v0.2.0', '5ebc990': 'v0.4.0', '5ed4b79': 'v0.6.0', '5ed986a': 'v0.6.0', '5f25624': 'v0.4.0',
    '5f3d292': 'v0.2.0', '5fd957b': 'v0.6.0', '60252e2': 'v0.6.0', '6046706': 'v0.4.0', '607342a': 'v0.2.0', '6098fea': 'v0.6.0',
    '609bb98': 'v0.6.0', '609fbcb': 'v0.6.0', '610614d': 'v0.5.0', '610c000': 'v0.6.0', '618c50f': 'v0.2.0', '6190914': 'v0.4.0',
    '61a0509': 'v0.6.0', '61d0464': 'v0.6.0', '6242e3b': 'v0.4.0', '628d640': 'v0.6.1', '62fe0b6': 'v0.4.0', '63078e8': 'v0.1.0',
    '63af4ed': 'v0.4.1', '63b3bf7': 'v0.6.0', '63d0a73': 'v0.6.0', '63f76da': 'v0.4.2', '64740cc': 'v0.6.0', '64fb78e': 'v0.4.0',
    '65b1491': 'v0.2.0', '65d940e': 'v0.3.0', '6611a7e': 'v0.1.0', '663e1a4': 'v0.4.0', '6663e4c': 'v0.1.0', '6679702': 'v0.6.0',
    '66a0ce1': 'v0.5.0', '66b0e27': 'v0.6.0', '66b3cdd': 'v0.6.1', '67352f1': 'v0.6.0', '6789f60': 'v0.1.0', '679c509': 'v0.4.0',
    '679c8ba': 'v0.6.0', '67c92c2': 'v0.6.0', '680a003': 'v0.3.0', '6836dac': 'v0.6.0', '6879af7': 'v0.3.0', '689cbec': 'v0.6.0',
    '68cd43d': 'v0.6.0', '68dab78': 'v0.1.0', '699f9a8': 'v0.1.0', '69bebb6': 'v0.6.0', '69d0938': 'v0.3.0', '69d82fd': 'v0.5.0',
    '69e9ffe': 'v0.4.0', '69f42b6': 'v0.4.0', '69f45a8': 'v0.6.0', '69f638f': 'v0.6.0', '6a74a0b': 'v0.6.0', '6aca0ca': 'v0.5.0',
    '6b09def': 'v0.6.0', '6b4f8c1': 'v0.6.0', '6c84f08': 'v0.4.0', '6cda4dc': 'v0.6.0', '6d38b9d': 'v0.5.0', '6d486dc': 'v0.6.1',
    '6df8bb8': 'v0.5.0', '6e52a11': 'v0.4.1', '6e9abb0': 'v0.6.1', '6ebe4d9': 'v0.6.0', '6ed1d38': 'v0.4.0', '6ed4510': 'v0.6.0',
    '6efa797': 'v0.6.0', '6f27a4c': 'v0.6.0', '700aebb': 'v0.5.0', '7027738': 'v0.6.0', '70eea6b': 'v0.5.0', '7175f4c': 'v0.6.0',
    '7182762': 'v0.3.0', '723513b': 'v0.4.0', '72cf4bb': 'v0.4.1', '732c39a': 'v0.6.0', '73de37c': 'v0.4.0', '745153f': 'v0.6.0',
    '752928e': 'v0.5.0', '75725e0': 'v0.6.0', '75bad96': 'v0.2.0', '75badde': 'v0.6.1', '75c9fa6': 'v0.4.0', '760a701': 'v0.3.0',
    '765821f': 'v0.6.0', '76a9d66': 'v0.6.0', '76bda19': 'v0.6.0', '76cdeef': 'v0.6.0', '76ce184': 'v0.6.1', '76ec298': 'v0.5.0',
    '77260e6': 'v0.6.0', '77956d3': 'v0.4.0', '779cff9': 'v0.1.0', '77b3046': 'v0.1.0', '77b3ba5': 'v0.6.0', '7831be6': 'v0.5.0',
    '78500ad': 'v0.6.0', '7864a63': 'v0.1.0', '78785eb': 'v0.1.0', '788644d': 'v0.6.0', '7889769': 'v0.6.0', '78ab26d': 'v0.6.0',
    '78b299a': 'v0.4.0', '78c36dd': 'v0.6.0', '78ce611': 'v0.6.0', '78cf1e6': 'v0.6.0', '791eb7b': 'v0.5.0', '7949263': 'v0.4.0',
    '79f375a': 'v0.6.0', '7a1177d': 'v0.6.0', '7aacec6': 'v0.1.0', '7acb6ab': 'v0.1.0', '7b58f48': 'v0.5.0', '7b76ff5': 'v0.6.1',
    '7bc831a': 'v0.6.0', '7beeb2a': 'v0.1.0', '7c29003': 'v0.6.0', '7c48d45': 'v0.5.0', '7c75bbd': 'v0.6.0', '7c81011': 'v0.4.1',
    '7c91356': 'v0.4.2', '7cfd531': 'v0.6.0', '7d33bbb': 'v0.3.0', '7d41473': 'v0.1.0', '7d53ec0': 'v0.1.0', '7d7c5b3': 'v0.4.1',
    '7da6a6a': 'v0.6.1', '7dd3732': 'v0.4.1', '7ddfa3b': 'v0.1.0', '7dea680': 'v0.6.0', '7df3a1e': 'v0.6.0', '7e070e7': 'v0.5.0',
    '7e22d51': 'v0.4.1', '7f033d6': 'v0.6.0', '7f39d3c': 'v0.6.0', '7f3b472': 'v0.1.0', '7f6b1d4': 'v0.6.0', '7f7eb56': 'v0.6.0',
    '7f9cc94': 'v0.6.0', '7faef0f': 'v0.1.0', '7fd1f23': 'v0.5.0', '8019263': 'v0.6.0', '80e9055': 'v0.6.0', '8135808': 'v0.5.0',
    '81971ac': 'v0.5.0', '81b1097': 'v0.6.0', '81b9aa0': 'v0.1.0', '81e4de7': 'v0.6.1', '8261590': 'v0.5.0', '826dd82': 'v0.4.0',
    '826f0c3': 'v0.1.0', '82e6c76': 'v0.4.0', '8303ebf': 'v0.1.0', '8345741': 'v0.6.0', '8381443': 'v0.4.2', '83c23e6': 'v0.6.0',
    '83ee019': 'v0.6.0', '8433368': 'v0.5.0', '845d493': 'v0.3.0', '84ba03e': 'v0.6.0', '84f5a52': 'v0.6.0', '855cba6': 'v0.3.0',
    '85aa3fe': 'v0.3.0', '85b9464': 'v0.6.0', '8608494': 'v0.3.0', '86314ef': 'v0.3.0', '864b283': 'v0.6.1', '86511ab': 'v0.5.0',
    '867b079': 'v0.2.0', '86bd642': 'v0.4.0', '86e000d': 'v0.3.0', '86ea127': 'v0.6.1', '8709494': 'v0.5.0', '874c322': 'v0.6.0',
    '875951b': 'v0.6.0', '879947b': 'v0.1.0', '87fb32d': 'v0.4.0', '8817c48': 'v0.4.0', '88dcc9c': 'v0.4.0', '8902190': 'v0.5.0',
    '8995d19': 'v0.1.0', '8a04060': 'v0.6.0', '8a898c1': 'v0.1.0', '8aea927': 'v0.6.0', '8af9cc8': 'v0.4.0', '8b433d0': 'v0.4.0',
    '8bbda47': 'v0.5.0', '8d26280': 'v0.1.0', '8d6546e': 'v0.6.0', '8d70f75': 'v0.6.0', '8d8519b': 'v0.2.0', '8d9b97b': 'v0.4.0',
    '8de0a9a': 'v0.6.0', '8df1faf': 'v0.6.0', '8e26e2d': 'v0.4.0', '8e7bdc9': 'v0.6.0', '8efd220': 'v0.6.0', '8f1b3ab': 'v0.4.0',
    '8f888d5': 'v0.4.2', '8fc779d': 'v0.4.0', '9000cf0': 'v0.5.0', '9032519': 'v0.6.0', '90a2a89': 'v0.5.0', '90bdb06': 'v0.4.0',
    '90cb303': 'v0.6.1', '914b321': 'v0.3.0', '9161793': 'v0.6.0', '91958f6': 'v0.3.0', '91ca465': 'v0.1.0', '9201f7b': 'v0.4.0',
    '92474fc': 'v0.4.1', '92a3081': 'v0.5.0', '92c77d5': 'v0.4.0', '92e2c67': 'v0.6.0', '9321ee0': 'v0.5.0', '936b4c1': 'v0.6.0',
    '9377f86': 'v0.5.0', '93dc883': 'v0.6.0', '940a0e5': 'v0.6.0', '945d792': 'v0.6.0', '949232c': 'v0.6.0', '94c4244': 'v0.2.0',
    '94ea3d7': 'v0.4.0', '94fc4c8': 'v0.6.1', '95509a7': 'v0.6.0', '9576fa0': 'v0.6.1', '957aa9f': 'v0.6.0', '9642682': 'v0.1.0',
    '96656e4': 'v0.1.0', '9682974': 'v0.5.0', '96a39e7': 'v0.1.0', '96a9600': 'v0.1.0', '970449b': 'v0.6.0', '9734c37': 'v0.2.0',
    '9760647': 'v0.1.0', '9770540': 'v0.3.0', '97b9921': 'v0.5.0', '97b9c92': 'v0.6.0', '9826dcf': 'v0.6.0', '98337c2': 'v0.3.0',
    '983e822': 'v0.2.0', '9849ff7': 'v0.5.0', '98abbb2': 'v0.4.0', '9994965': 'v0.6.0', '9a0f41c': 'v0.3.0', '9a1ca4c': 'v0.6.0',
    '9a818c1': 'v0.6.0', '9a93de5': 'v0.3.0', '9ae3215': 'v0.6.0', '9af8d9d': 'v0.1.0', '9b45ec9': 'v0.1.0', '9be2bd3': 'v0.6.0',
    '9c9e87e': 'v0.4.1', '9ccf5cd': 'v0.4.0', '9cd2764': 'v0.1.0', '9d14971': 'v0.6.1', '9d4dda6': 'v0.3.0', '9d7bdd8': 'v0.6.0',
    '9db3122': 'v0.6.0', '9ddc8a4': 'v0.4.2', '9df6351': 'v0.4.0', '9df6ba0': 'v0.1.0', '9e255df': 'v0.3.0', '9e2c0c4': 'v0.1.0',
    '9e902f8': 'v0.4.0', '9eb180b': 'v0.1.0', '9ece053': 'v0.6.0', '9ee032f': 'v0.4.0', '9efe6a4': 'v0.1.0', '9f09f5c': 'v0.2.0',
    '9f3531b': 'v0.4.0', '9f4e8b1': 'v0.2.0', '9f76718': 'v0.1.0', '9f7a382': 'v0.4.0', '9f7b1d8': 'v0.4.2', '9f8c49e': 'v0.6.0',
    'a00f310': 'v0.3.0', 'a02cc42': 'v0.6.0', 'a04c234': 'v0.5.0', 'a1017f1': 'v0.6.0', 'a1418c7': 'v0.6.0', 'a17aa68': 'v0.4.0',
    'a19f759': 'v0.1.0', 'a1a8fd3': 'v0.4.0', 'a1c9900': 'v0.5.0', 'a1ddaf0': 'v0.5.0', 'a1f3a2c': 'v0.6.1', 'a26fa5b': 'v0.3.0',
    'a2b4eca': 'v0.6.1', 'a2be9e8': 'v0.6.0', 'a40b5ad': 'v0.6.0', 'a46c43e': 'v0.6.0', 'a4df53c': 'v0.6.0', 'a4e4e22': 'v0.6.0',
    'a54982f': 'v0.6.1', 'a551217': 'v0.6.0', 'a6428e0': 'v0.4.1', 'a6b698d': 'v0.5.0', 'a6dabd7': 'v0.6.0', 'a71965b': 'v0.1.0',
    'a753673': 'v0.5.0', 'a76d386': 'v0.1.0', 'a786ca9': 'v0.4.2', 'a7b607d': 'v0.6.0', 'a84798c': 'v0.6.1', 'a86b258': 'v0.6.1',
    'a86bd33': 'v0.1.0', 'a892689': 'v0.4.0', 'a8c7917': 'v0.1.0', 'a8cc50f': 'v0.1.0', 'a907aa3': 'v0.4.2', 'a91d514': 'v0.4.0',
    'a925c6d': 'v0.3.0', 'a9305dc': 'v0.6.0', 'a98dd0f': 'v0.4.1', 'a9ed0d5': 'v0.5.0', 'a9f881c': 'v0.6.0', 'aa0eeb2': 'v0.4.0',
    'aa1d59a': 'v0.6.0', 'aa2d495': 'v0.5.0', 'aa8403f': 'v0.4.0', 'aadde47': 'v0.6.1', 'aafd571': 'v0.3.0', 'ab72bad': 'v0.2.0',
    'ab7507c': 'v0.6.0', 'ab7a21a': 'v0.1.0', 'ab80a9b': 'v0.1.0', 'ab9cc50': 'v0.6.0', 'abc9679': 'v0.6.0', 'abeac58': 'v0.4.0',
    'abff1b2': 'v0.5.0', 'ac27843': 'v0.6.0', 'ac69f28': 'v0.6.0', 'ac77229': 'v0.6.0', 'ac801c3': 'v0.5.0', 'ad10be4': 'v0.6.1',
    'ad1870a': 'v0.3.0', 'ad2c3d5': 'v0.2.0', 'ad39db8': 'v0.6.1', 'ad48d96': 'v0.4.0', 'ad6af40': 'v0.5.0', 'ad85700': 'v0.6.0',
    'adf3908': 'v0.4.0', 'ae27b4f': 'v0.1.0', 'ae3214f': 'v0.4.0', 'ae62bae': 'v0.6.0', 'ae8965d': 'v0.6.1', 'aea5acd': 'v0.6.0',
    'aed5afc': 'v0.1.0', 'aef87b2': 'v0.4.0', 'af32cd4': 'v0.6.1', 'af69e68': 'v0.3.0', 'af74824': 'v0.6.0', 'af78ca6': 'v0.6.1',
    'afc09ad': 'v0.6.0', 'afcf7c2': 'v0.6.0', 'b030fb3': 'v0.4.0', 'b082dbc': 'v0.3.0', 'b0a6acd': 'v0.6.0', 'b0e8b18': 'v0.3.0',
    'b108e3d': 'v0.3.0', 'b125251': 'v0.1.0', 'b16aeaa': 'v0.6.0', 'b17498c': 'v0.6.0', 'b1dcb3c': 'v0.4.2', 'b268bb1': 'v0.6.0',
    'b331f80': 'v0.2.0', 'b399537': 'v0.6.0', 'b4029bc': 'v0.3.0', 'b432b8a': 'v0.4.0', 'b440752': 'v0.6.1', 'b47f8d8': 'v0.6.0',
    'b480266': 'v0.6.0', 'b4a3a4a': 'v0.6.0', 'b53ddd0': 'v0.5.0', 'b558560': 'v0.6.0', 'b56e032': 'v0.5.0', 'b58f4d5': 'v0.6.0',
    'b60a3dd': 'v0.6.1', 'b60db29': 'v0.6.0', 'b625c36': 'v0.4.2', 'b63eef9': 'v0.6.0', 'b677cd0': 'v0.6.0', 'b67d51c': 'v0.6.0',
    'b68971c': 'v0.3.0', 'b6924f5': 'v0.6.0', 'b69fce3': 'v0.2.0', 'b7073a0': 'v0.4.0', 'b7c8ba2': 'v0.5.0', 'b7cce9f': 'v0.6.0',
    'b7f5cb2': 'v0.6.0', 'b7f7acf': 'v0.6.0', 'b835218': 'v0.1.0', 'b88c399': 'v0.6.0', 'b8bb150': 'v0.6.0', 'b903d2c': 'v0.1.0',
    'b930653': 'v0.6.0', 'b9327b0': 'v0.4.0', 'b940014': 'v0.4.1', 'b9559e4': 'v0.1.0', 'ba2e849': 'v0.6.1', 'ba6d7b0': 'v0.3.0',
    'bac5ace': 'v0.5.0', 'bad6c33': 'v0.6.1', 'bae2b3d': 'v0.5.0', 'bafe3a3': 'v0.3.0', 'bb174a4': 'v0.6.0', 'bb385ba': 'v0.5.0',
    'bb81618': 'v0.1.0', 'bb8f3ce': 'v0.6.0', 'bb9964c': 'v0.6.0', 'bbc232a': 'v0.6.0', 'bbc44c6': 'v0.3.0', 'bbca07d': 'v0.2.0',
    'bc67583': 'v0.5.0', 'bc7cb33': 'v0.6.0', 'bc95160': 'v0.6.0', 'bca09cf': 'v0.4.0', 'bcc5955': 'v0.6.0', 'bcedff8': 'v0.3.0',
    'bd3710a': 'v0.6.1', 'bd57603': 'v0.2.0', 'bda9f08': 'v0.1.0', 'bdf0ff9': 'v0.3.0', 'be11b30': 'v0.6.0', 'be7b271': 'v0.6.1',
    'bebb228': 'v0.3.0', 'bebbab8': 'v0.5.0', 'bedf959': 'v0.6.0', 'bf45458': 'v0.6.0', 'bf84ba5': 'v0.6.1', 'bfad3b2': 'v0.3.0',
    'c01ab9d': 'v0.6.0', 'c02dfe9': 'v0.4.0', 'c02f3c6': 'v0.4.0', 'c03ce8b': 'v0.2.0', 'c045c98': 'v0.6.0', 'c0dcf17': 'v0.6.0',
    'c112723': 'v0.2.0', 'c1ddc1d': 'v0.6.0', 'c2043d9': 'v0.6.0', 'c311f88': 'v0.1.0', 'c405c15': 'v0.1.0', 'c45123a': 'v0.6.0',
    'c4645de': 'v0.6.0', 'c5497ed': 'v0.6.0', 'c55e7f6': 'v0.6.1', 'c5810ad': 'v0.2.0', 'c583c59': 'v0.6.0', 'c5a4a3d': 'v0.2.0',
    'c614de4': 'v0.5.0', 'c61a3e4': 'v0.4.0', 'c61e10e': 'v0.4.0', 'c625e76': 'v0.1.0', 'c6a3a5d': 'v0.6.0', 'c6bdd00': 'v0.4.0',
    'c759fd2': 'v0.4.1', 'c773638': 'v0.3.0', 'c7c6a8f': 'v0.1.0', 'c7e159d': 'v0.6.0', 'c7e8f1c': 'v0.6.0', 'c7f403b': 'v0.4.0',
    'c83113a': 'v0.2.0', 'c855b69': 'v0.6.0', 'c89022b': 'v0.4.0', 'c8e4761': 'v0.4.0', 'c8ea6db': 'v0.6.0', 'c925b58': 'v0.4.0',
    'c926c7b': 'v0.1.0', 'c93075e': 'v0.4.2', 'c964446': 'v0.4.0', 'c966712': 'v0.6.0', 'c9ad741': 'v0.2.0', 'c9ca060': 'v0.4.0',
    'c9d44ad': 'v0.6.0', 'ca07244': 'v0.4.1', 'ca0d60c': 'v0.2.0', 'ca16b85': 'v0.4.0', 'ca3b5b4': 'v0.5.0', 'ca41bd7': 'v0.4.0',
    'ca4cdda': 'v0.2.0', 'cacc425': 'v0.2.0', 'cae18dc': 'v0.5.0', 'cae3ff7': 'v0.6.0', 'caf3979': 'v0.6.0', 'caff615': 'v0.5.0',
    'cb43ada': 'v0.5.0', 'cb786d0': 'v0.1.0', 'cbb1987': 'v0.6.0', 'cc08e67': 'v0.3.0', 'cc447e1': 'v0.6.0', 'cd2e84a': 'v0.6.0',
    'cd57df3': 'v0.2.0', 'cd7913a': 'v0.4.0', 'cd99c92': 'v0.6.0', 'cdcfda4': 'v0.6.0', 'ce18514': 'v0.6.0', 'ce2b7b0': 'v0.6.0',
    'cf20f5b': 'v0.1.0', 'cf4334d': 'v0.3.0', 'cf5931a': 'v0.6.0', 'cf8d4a2': 'v0.6.0', 'cf9e7d8': 'v0.2.0', 'cfacaf1': 'v0.5.0',
    'cfdc1f4': 'v0.5.0', 'cff3626': 'v0.4.0', 'd0170f6': 'v0.4.0', 'd0859f1': 'v0.1.0', 'd0aea39': 'v0.4.0', 'd12d830': 'v0.6.0',
    'd13e499': 'v0.6.1', 'd15faf9': 'v0.6.0', 'd187f71': 'v0.6.0', 'd18e838': 'v0.2.0', 'd1a7fc1': 'v0.6.0', 'd1afebb': 'v0.5.0',
    'd276b22': 'v0.6.0', 'd2af40b': 'v0.6.0', 'd2c033a': 'v0.3.0', 'd2c85ba': 'v0.1.0', 'd30ab8a': 'v0.5.0', 'd312257': 'v0.6.0',
    'd3886f7': 'v0.3.0', 'd3944ef': 'v0.6.1', 'd40623d': 'v0.3.0', 'd41f7fc': 'v0.5.0', 'd445bb0': 'v0.5.0', 'd48690a': 'v0.6.0',
    'd4d17ab': 'v0.5.0', 'd53067a': 'v0.6.0', 'd534f1b': 'v0.4.0', 'd53586f': 'v0.6.1', 'd53d20c': 'v0.1.0', 'd5718bd': 'v0.6.0',
    'd5ae801': 'v0.5.0', 'd5f497e': 'v0.6.0', 'd628038': 'v0.6.0', 'd6a581a': 'v0.4.0', 'd6d15a0': 'v0.1.0', 'd6d8719': 'v0.3.0',
    'd70e8b0': 'v0.6.0', 'd741f7b': 'v0.3.0', 'd74418d': 'v0.1.0', 'd76a82d': 'v0.6.0', 'd7c80bd': 'v0.5.0', 'd8101ff': 'v0.2.0',
    'd810fb2': 'v0.6.0', 'd83e970': 'v0.3.0', 'd94aa0c': 'v0.3.0', 'da3fab3': 'v0.1.0', 'da53daa': 'v0.6.0', 'da99bba': 'v0.6.0',
    'dacff69': 'v0.4.0', 'db0535a': 'v0.4.0', 'db271c6': 'v0.1.0', 'db72554': 'v0.6.0', 'db73835': 'v0.1.0', 'db9b219': 'v0.6.0',
    'dba0dcc': 'v0.3.0', 'dc0bffe': 'v0.6.0', 'dc19182': 'v0.6.0', 'dc20a59': 'v0.6.1', 'dc54ca6': 'v0.6.1', 'dc5abf8': 'v0.5.0',
    'dc70a4a': 'v0.6.0', 'dce0f84': 'v0.1.0', 'dce7921': 'v0.5.0', 'dddd9aa': 'v0.6.1', 'de1ab58': 'v0.4.2', 'de9c94f': 'v0.1.0',
    'deb6273': 'v0.5.0', 'df33b3e': 'v0.5.0', 'df381c5': 'v0.4.0', 'df3a2ee': 'v0.4.0', 'df40c1b': 'v0.2.0', 'df643df': 'v0.2.0',
    'df871e1': 'v0.1.0', 'df88dec': 'v0.4.0', 'df9bd75': 'v0.1.0', 'dfaf6f6': 'v0.6.0', 'dfc14c0': 'v0.3.0', 'dfdf823': 'v0.5.0',
    'dfdfd18': 'v0.1.0', 'e05aff6': 'v0.4.0', 'e088713': 'v0.5.0', 'e0d4489': 'v0.4.0', 'e0e9f8a': 'v0.1.0', 'e0f0ac5': 'v0.6.0',
    'e13f284': 'v0.6.0', 'e145f96': 'v0.6.0', 'e16967b': 'v0.4.0', 'e17dd01': 'v0.6.0', 'e1a3ff7': 'v0.1.0', 'e1aa3bb': 'v0.4.0',
    'e24124d': 'v0.4.0', 'e250a12': 'v0.6.0', 'e269fcf': 'v0.2.0', 'e26d7a1': 'v0.5.0', 'e2d6523': 'v0.5.0', 'e2d87d3': 'v0.6.0',
    'e2d91ed': 'v0.6.0', 'e2fb7b6': 'v0.4.2', 'e318dce': 'v0.1.0', 'e388e59': 'v0.6.0', 'e3a63a3': 'v0.4.0', 'e42f360': 'v0.6.1',
    'e43e40e': 'v0.5.0', 'e47c819': 'v0.1.0', 'e4a6024': 'v0.6.0', 'e4bb0fe': 'v0.5.0', 'e4f5498': 'v0.4.0', 'e53669d': 'v0.6.0',
    'e573604': 'v0.4.0', 'e59b47d': 'v0.4.0', 'e5aa59f': 'v0.6.0', 'e66de12': 'v0.6.0', 'e6bd1b2': 'v0.6.1', 'e6f2b06': 'v0.6.0',
    'e7e38fb': 'v0.4.0', 'e80ea6c': 'v0.1.0', 'e90871c': 'v0.1.0', 'e91abb0': 'v0.6.0', 'e928aeb': 'v0.1.0', 'e92d4f1': 'v0.4.0',
    'e94d163': 'v0.4.0', 'e9875fc': 'v0.6.0', 'e99d9b7': 'v0.6.0', 'e9b1912': 'v0.3.0', 'e9b3076': 'v0.5.0', 'e9cdc3c': 'v0.6.0',
    'e9d2230': 'v0.3.0', 'ea00134': 'v0.4.2', 'ea492c7': 'v0.6.0', 'eb03105': 'v0.6.0', 'eb29dfc': 'v0.4.0', 'eb41765': 'v0.6.0',
    'eb48817': 'v0.3.0', 'eb5a151': 'v0.6.0', 'eb78bb7': 'v0.4.2', 'eb7e7c1': 'v0.6.0', 'ec041c9': 'v0.2.0', 'ec14e7a': 'v0.6.0',
    'ec64479': 'v0.6.0', 'ec6be60': 'v0.3.0', 'ed16a31': 'v0.5.0', 'ed2e7fd': 'v0.6.0', 'ed4a981': 'v0.6.1', 'ed4b0be': 'v0.4.0',
    'ed6a469': 'v0.6.0', 'ed713e4': 'v0.6.0', 'ed74a93': 'v0.6.0', 'edae95d': 'v0.1.0', 'edb3e7f': 'v0.6.0', 'edf6f46': 'v0.1.0',
    'ee06f3d': 'v0.6.0', 'ee77c88': 'v0.6.0', 'ef12b6a': 'v0.6.0', 'ef482df': 'v0.6.1', 'efc9861': 'v0.6.0', 'eff6039': 'v0.6.1',
    'f008edb': 'v0.6.0', 'f013f61': 'v0.1.0', 'f036139': 'v0.6.0', 'f0f6402': 'v0.1.0', 'f113dc1': 'v0.6.0', 'f11624b': 'v0.6.0',
    'f14f5e0': 'v0.3.0', 'f16b74b': 'v0.6.0', 'f18cae6': 'v0.4.0', 'f1e0035': 'v0.6.0', 'f1ef634': 'v0.1.0', 'f22230d': 'v0.4.0',
    'f234896': 'v0.4.1', 'f245e4f': 'v0.6.1', 'f2470f2': 'v0.1.0', 'f27beae': 'v0.3.0', 'f390ff6': 'v0.1.0', 'f3e569f': 'v0.3.0',
    'f3f9639': 'v0.3.0', 'f455385': 'v0.4.0', 'f45b9a2': 'v0.6.0', 'f468f26': 'v0.6.0', 'f4a6c43': 'v0.1.0', 'f4c387b': 'v0.6.0',
    'f51065a': 'v0.5.0', 'f54c60e': 'v0.4.0', 'f553aec': 'v0.5.0', 'f57facc': 'v0.6.0', 'f5873ff': 'v0.1.0', 'f5bfc4c': 'v0.4.2',
    'f6022b7': 'v0.4.0', 'f62243b': 'v0.4.0', 'f6a1e17': 'v0.5.0', 'f6ae14c': 'v0.6.1', 'f6dfbfa': 'v0.6.0', 'f726c49': 'v0.2.0',
    'f749590': 'v0.4.0', 'f795803': 'v0.1.0', 'f7990ec': 'v0.5.0', 'f82b605': 'v0.3.0', 'f88986a': 'v0.1.0', 'f8d7286': 'v0.3.0',
    'f90bc2f': 'v0.6.0', 'f9147f4': 'v0.6.0', 'f936658': 'v0.1.0', 'f937c9a': 'v0.1.0', 'f941e11': 'v0.4.2', 'f94c040': 'v0.3.0',
    'f95d50b': 'v0.6.0', 'f9d3734': 'v0.2.0', 'fa45f1b': 'v0.6.0', 'fa4cd8e': 'v0.6.1', 'fa4fe45': 'v0.6.0', 'fa71b41': 'v0.6.0',
    'fa755d5': 'v0.1.0', 'fa8fae5': 'v0.6.1', 'faa1d4e': 'v0.6.0', 'fac0550': 'v0.5.0', 'fb0a377': 'v0.4.0', 'fb4690b': 'v0.5.0',
    'fb4b2fc': 'v0.5.0', 'fb61199': 'v0.6.0', 'fb8bfb1': 'v0.4.0', 'fbcb540': 'v0.5.0', 'fbdb27f': 'v0.6.0', 'fc33b0b': 'v0.1.0',
    'fc7dfa6': 'v0.3.0', 'fc9e3b8': 'v0.3.0', 'fca32a7': 'v0.1.0', 'fcb2924': 'v0.6.0', 'fd37b2e': 'v0.3.0', 'fd55568': 'v0.6.0',
    'fd688fd': 'v0.6.0', 'fd82fb4': 'v0.3.0', 'fd8486b': 'v0.6.0', 'fe1e2d6': 'v0.5.0', 'fed43e5': 'v0.6.0',
  };

  function resolveTagFamily(row) {
    const ver = row.prestus_ver || '';
    // Full describe format: v1.2.3-N-gHASH → strip the commit suffix
    const describeMatch = ver.match(/^(v[\d.]+(?:-[a-zA-Z]\w*)*)-\d+-g[0-9a-f]+$/);
    if (describeMatch) return describeMatch[1];
    // Clean tag: v1.2.3 or v1.2.3-rc1
    if (ver.match(/^v[\d.]+(?:-[a-zA-Z]\w*)*$/)) return ver;
    // Legacy: prestus_ver was a bare hash, or fall back to prestus_hash
    const legacyHash = (row.prestus_hash || ver || '').slice(0, 7);
    return HASH_INFO[legacyHash] || 'unknown';
  }

  function resolveVersion(row) {
    const hash = (row.prestus_hash || '').slice(0, 7);
    const ver  = row.prestus_ver  || '';

    // Full git-describe format: v1.2.3-N-gABCDEF (non-tag commit)
    const describeMatch = ver.match(/^(v[\d.]+(?:-\w+)*)-(\d+)-g([0-9a-f]{7,})$/);
    if (describeMatch) {
      const [, tag, n, h] = describeMatch;
      return `${tag} - ${n} - ${h.slice(0, 7)}`;
    }

    // Clean tag commit: prestus_ver is just the tag (e.g. "v0.6.1")
    if (ver.match(/^v[\d.]+(?:-\w+)*$/)) {
      return hash ? `${ver} - 0 - ${hash}` : ver;
    }

    // Legacy: look up by hash
    const legacyHash = (hash || ver || '').slice(0, 7);
    if (legacyHash && HASH_INFO[legacyHash]) {
      return `${HASH_INFO[legacyHash]} - ${legacyHash}`;
    }

    return legacyHash || 'unknown';
  }

  function classifyMedium(row) {
    const med = row.medium || 'unknown';
    if (med === 'layered' && row.pct_enabled) return 'layered + pCT';
    return med;
  }

  function fg()   { return getComputedStyle(document.body).getPropertyValue('--md-default-fg-color')         || '#333'; }
  function fgL()  { return getComputedStyle(document.body).getPropertyValue('--md-default-fg-color--light')  || '#999'; }
  function gridC(){ return getComputedStyle(document.body).getPropertyValue('--md-default-fg-color--lightest')|| '#e0e0e0'; }

  function countBy(arr, fn) {
    return arr.reduce((a, r) => { const v = fn(r); a[v] = (a[v]||0)+1; return a; }, {});
  }

  function doughnutChart(id, labels, data) {
    new Chart(document.getElementById(id), {
      type: 'doughnut',
      data: { labels, datasets: [{ data, backgroundColor: PAL.slice(0, labels.length), borderWidth: 2 }] },
      options: { plugins: { legend: { position: 'right', labels: { color: fg(), boxWidth: 12, font: { size: 11 } } } }, cutout: '60%' }
    });
  }

  async function init() {
    const client = supabase.createClient(URL, KEY);

    const { data, error } = await client
      .from('events')
      .select('received_at,event,status,duration_s,uuid,sim_platform,code_type,pipeline_mode,medium,pct_enabled,prestus_ver,prestus_hash')
      .in('event', ['run_end', 'run_error'])
      .limit(5000);

    document.getElementById('telem-loading').style.display = 'none';

    if (error) {
      const el = document.getElementById('telem-error');
      el.textContent = 'Could not load telemetry data: ' + error.message;
      el.style.display = 'block';
      return;
    }

    // ── Stat cards ─────────────────────────────────────────────────
    const total     = data.length;
    const layeredRuns = data.filter(r => classifyMedium(r).startsWith('layered'));
    const layeredSucc = layeredRuns.filter(r => r.status === 'success');
    const successPct  = layeredRuns.length ? Math.round(layeredSucc.length / layeredRuns.length * 100) + '%' : '—';
    const users     = new Set(data.map(r => r.uuid).filter(Boolean)).size;
    const versions  = new Set(data.map(resolveVersion).filter(v => v !== 'unknown')).size;

    const statsEl = document.getElementById('telem-stats');
    statsEl.style.display = 'grid';
    [
      ['Completed runs', total],
      ['Success rate (layered)', successPct],
      ['Versions',       versions],
      ['Unique users',   users],
    ].forEach(([label, val]) => {
      statsEl.insertAdjacentHTML('beforeend',
        `<div class="telem-stat"><div class="tv">${val}</div><div class="tl">${label}</div></div>`);
    });

    document.getElementById('telem-charts').style.display = 'grid';

    // ── Simulations per month, stacked by tag family ───────────────
    const allVersions = [...new Set(data.map(resolveTagFamily))].filter(v => v !== 'unknown').sort((a, b) => a === 'untagged' ? -1 : b === 'untagged' ? 1 : a.localeCompare(b));
    const monthSet = new Set(data.map(r => r.received_at?.slice(0, 7)).filter(Boolean));
    const months   = [...monthSet].sort();

    const mvCount = {};
    data.forEach(r => {
      const key = r.received_at?.slice(0, 7) + '|' + resolveTagFamily(r);
      if (key) mvCount[key] = (mvCount[key] || 0) + 1;
    });

    new Chart(document.getElementById('tc-timeline'), {
      type: 'bar',
      data: {
        labels: months,
        datasets: allVersions.map((ver, i) => ({
          label: ver,
          data: months.map(m => mvCount[m + '|' + ver] || 0),
          backgroundColor: PAL[i % PAL.length],
          borderRadius: 2,
          borderSkipped: false,
        }))
      },
      options: {
        plugins: {
          legend: {
            position: 'right',
            maxWidth: 160,
            labels: {
              color: fg(), boxWidth: 10, font: { size: 10 },
              generateLabels: chart => Chart.defaults.plugins.legend.labels.generateLabels(chart)
                .map(l => ({ ...l, text: l.text.length > 18 ? l.text.slice(0, 18) + '…' : l.text }))
            }
          }
        },
        scales: {
          x: { stacked: true, ticks: { color: fgL(), maxRotation: 45, font: { size: 10 } }, grid: { display: false } },
          y: { stacked: true, ticks: { color: fgL() }, grid: { color: gridC() }, beginAtZero: true }
        }
      }
    });

    // ── Execution platform ─────────────────────────────────────────
    const pc = countBy(data, r => r.sim_platform ?? '(unknown)');
    doughnutChart('tc-platform', Object.keys(pc), Object.values(pc));

    // ── Execution backend (code_type) ──────────────────────────────
    const cc = countBy(data, r => r.code_type ?? '(unknown)');
    const ck = ['matlab_cpu','matlab_gpu','cpp_cpu','cpp_gpu'].filter(k => cc[k]);
    const unknownCC = Object.keys(cc).filter(k => !['matlab_cpu','matlab_gpu','cpp_cpu','cpp_gpu'].includes(k));
    doughnutChart('tc-codetype',
      [...ck, ...unknownCC],
      [...ck.map(k => cc[k]), ...unknownCC.map(k => cc[k])]
    );

    // ── Pipeline mode ──────────────────────────────────────────────
    const pmc = countBy(data, r => r.pipeline_mode ?? '(unknown)');
    doughnutChart('tc-pipemode', Object.keys(pmc), Object.values(pmc));

    // ── Simulation medium (layered+pCT as its own category) ────────
    const mc = countBy(data, classifyMedium);
    doughnutChart('tc-medium', Object.keys(mc), Object.values(mc));

    // ── Mean run duration by medium (all completed runs) ─────────
    const durByMedium = {};
    data.filter(r => r.duration_s).forEach(r => {
      const med = classifyMedium(r);
      if (!durByMedium[med]) durByMedium[med] = [];
      durByMedium[med].push(r.duration_s / 60);
    });
    const durOrder  = ['water', 'phantom', 'layered', 'layered + pCT'];
    const durLabels = [...durOrder.filter(m => durByMedium[m]),
                       ...Object.keys(durByMedium).filter(m => !durOrder.includes(m))];
    const durMeans  = durLabels.map(m => {
      const vals = durByMedium[m];
      return +(vals.reduce((a, b) => a + b, 0) / vals.length).toFixed(1);
    });
    new Chart(document.getElementById('tc-duration'), {
      type: 'bar',
      data: {
        labels: durLabels,
        datasets: [{
          data: durMeans,
          backgroundColor: durLabels.map((_, i) => PAL[i % PAL.length]),
          borderRadius: 4,
          borderSkipped: false,
        }]
      },
      options: {
        plugins: { legend: { display: false } },
        scales: {
          x: { ticks: { color: fgL(), font: { size: 11 } }, grid: { display: false } },
          y: { ticks: { color: fgL(), callback: v => v + ' min' }, grid: { color: gridC() }, beginAtZero: true }
        }
      }
    });
    // ── Success rate by medium ─────────────────────────────────────
    const srByMedium = {};
    data.forEach(r => {
      const med = classifyMedium(r);
      if (!srByMedium[med]) srByMedium[med] = { ok: 0, total: 0 };
      srByMedium[med].total++;
      if (r.status === 'success') srByMedium[med].ok++;
    });
    const srOrder  = ['water', 'phantom', 'layered', 'layered + pCT'];
    const srLabels = [...srOrder.filter(m => srByMedium[m]),
                      ...Object.keys(srByMedium).filter(m => !srOrder.includes(m))];
    const srValues = srLabels.map(m => +(srByMedium[m].ok / srByMedium[m].total * 100).toFixed(1));
    new Chart(document.getElementById('tc-success-rate'), {
      type: 'bar',
      data: {
        labels: srLabels,
        datasets: [{
          data: srValues,
          backgroundColor: srLabels.map((_, i) => PAL[i % PAL.length]),
          borderRadius: 4,
          borderSkipped: false,
        }]
      },
      options: {
        plugins: { legend: { display: false } },
        scales: {
          x: { ticks: { color: fgL(), font: { size: 11 } }, grid: { display: false } },
          y: { min: 0, max: 100, ticks: { color: fgL(), callback: v => v + '%' }, grid: { color: gridC() } }
        }
      }
    });
  }

  if (typeof Chart !== 'undefined' && typeof supabase !== 'undefined') {
    init();
  } else {
    window.addEventListener('load', init);
  }
})();
</script>
