function drawSingleOP(data) {
  // let oldPlot = document.querySelector("#plot");
  // if (oldPlot) {
  //   oldPlot.remove();
  // }
  let fig = document.createElement("div");
  fig.setAttribute("id", `fig${window.app.resultIndices.slice(-1)[0]}`);
  fig.setAttribute("style", "margin-left: auto; margin-right: auto;");
  document.querySelector("#chart").appendChild(fig);

  if (!(typeof data === "number")) {
    let t = data.t;
    let ir = data.ir;
    let vr = data.vr;
    let im = data.im;
    let vhb = data.vhb;

    const padding = 30;
    const width = 400;
    const height = 300;

    const svg = d3
      .select(`#fig${window.app.resultIndices.slice(-1)[0]}`)
      .append("svg")
      .attr("id", `plot${window.app.resultIndices.slice(-1)[0]}`)
      .attr("width", width)
      .attr("height", height)
      .attr("style", "margin-left: auto; margin-right: auto;");
    const tScale = d3
      .scaleLinear()
      .domain([0, Math.max(...t) / 1e-6])
      .range([padding, width - padding]);
    const iScale = d3
      .scaleLinear()
      .domain([Math.min(...ir), Math.max(...ir)])
      .range([height - padding, padding]);
    const vScale = d3
      .scaleLinear()
      .domain([Math.min(...vr, ...vhb), Math.max(...vr, ...vhb)])
      .range([height - padding, padding]);
    const currentLine = d3
      .line()
      .x((d) => tScale(d.x))
      .y((d) => iScale(d.y))
      .curve(d3.curveLinear);
    const voltageLine = d3
      .line()
      .x((d) => tScale(d.x))
      .y((d) => vScale(d.y))
      .curve(d3.curveLinear);

    let irData = t.map((e, k) => ({ x: e / 1e-6, y: ir[k] }));
    let imData = t.map((e, k) => ({ x: e / 1e-6, y: im[k] }));
    let vrData = t.map((e, k) => ({ x: e / 1e-6, y: vr[k] }));
    let vhbData = t.map((e, k) => ({ x: e / 1e-6, y: vhb[k] }));

    svg
      .append("path")
      .classed("data", true)
      .attr("d", currentLine(irData))
      .style("fill", "none")
      .style("stroke", "blue");
    svg
      .append("path")
      .classed("data", true)
      .attr("d", currentLine(imData))
      .style("fill", "none")
      .style("stroke", "red");
    svg
      .append("path")
      .classed("data", true)
      .attr("d", voltageLine(vrData))
      .style("fill", "none")
      .style("stroke", "orange");
    svg
      .append("path")
      .classed("data", true)
      .attr("d", voltageLine(vhbData))
      .style("fill", "none")
      .style("stroke", "darkgreen");

    // svg.call(
    //   d3.zoom().on("zoom", () => {
    //     svg.attr(
    //       "transform",
    //       `translate(${d3.event.translate}) scale(${d3.event.scale})`
    //     );
    //   })
    // );

    let tAxis = d3.axisBottom(tScale);
    let iAxis = d3.axisLeft(iScale);
    let vAxis = d3.axisRight(vScale);
    svg
      .append("g")
      .attr("class", "t-axis")
      .attr("transform", `translate(0, ${height - padding})`)
      .call(tAxis);
    svg
      .selectAll("g.t-axis g.tick")
      .append("line")
      .classed("grid-line", true)
      .attr("y2", -(height - 2 * padding))
      .attr("stroke", "#dddddd");
    svg.append("g").attr("class", "i-axis").attr("transform", `translate(${padding})`).call(iAxis);
    svg
      .selectAll("g.i-axis g.tick")
      .append("line")
      .classed("grid-line", true)
      .attr("x2", width - 2 * padding)
      .attr("stroke", "#dddddd");
    svg
      .append("g")
      .attr("class", "v-axis")
      .attr("transform", `translate(${width - padding})`)
      .call(vAxis);
    svg
      .selectAll("g.v-axis g.tick")
      .append("line")
      .classed("grid-line", true)
      .attr("x2", width - 2 * padding)
      .attr("stroke", "#dddddd");
  } else {
    let errorValue = data;
    let errorMessage = document.createElement("p");
    errorMessage.classList.add("result-error-message");
    errorMessage.innerHTML = `Steady-state could not be found. The maximum possible output power is ${errorValue}.`;
    fig.appendChild(errorMessage);
  }

  let caption = document.createElement("p");
  caption.innerHTML = `Figure ${window.app.resultIndices.slice(-1)[0]}`;
  caption.setAttribute("style", "text-align: center;");
  fig.appendChild(caption);
}

function drawDesign(x, ys, parameters) {
  for (const code of ["00", "01", "10", "11"]) {
    let old = document.querySelector(`#fig${code}`);
    if (old) {
      old.remove();
    }
  }

  let y00;
  let y01;
  let y10;
  let y11;
  [y00, y01, y10, y11] = ys;

  let fig00 = document.createElement("div");
  fig00.setAttribute("id", "fig00");
  fig00.setAttribute("style", "margin-left: auto; margin-right:auto;");
  let fig01 = document.createElement("div");
  fig01.setAttribute("id", "fig01");
  fig01.setAttribute("style", "margin-left: auto; margin-right:auto;");
  let fig10 = document.createElement("div");
  fig10.setAttribute("id", "fig10");
  fig10.setAttribute("style", "margin-left: auto; margin-right:auto;");
  let fig11 = document.createElement("div");
  fig11.setAttribute("id", "fig11");
  fig11.setAttribute("style", "margin-left: auto; margin-right:auto;");
  document.querySelector("#chart").appendChild(fig00);
  document.querySelector("#chart").appendChild(fig01);
  document.querySelector("#chart").appendChild(fig10);
  document.querySelector("#chart").appendChild(fig11);

  drawDesignChart1(x, y00, "00", "$I_{\\rm rms,pri}$ / A", "Primary RMS Current", 1e-6, 1);
  drawDesignChart1(x, y01, "01", "$\\Psi_{\\rm pri}$ / mWb", "Magnetic Flux linked to Primary", 1e-6, 1e-3);
  drawDesignChart1(x, y10, "10", "$Q_{\\rm comm}$ / nC", "Electric Charge for Low-High Commutation", 1e-6, 1e-9);
  drawDesignChart1(x, y11, "11", "$f_{\\rm sw}$ / kHz", "Switching Frequency", 1e-6, 1e3);
}

function drawDesignChart1(x, y, pos, ylabel = "", caption = "", kx = 1e-6, ky = 1) {
  const leftPadding = 60;
  const bottomPadding = 40;
  const rightPadding = 10;
  const topPadding = 10;
  const width = 400;
  const height = 300;
  const colors = ["red", "green", "blue"];

  const svg = d3
    .select(`#fig${pos}`)
    .append("svg")
    .attr("width", width)
    .attr("height", height)
    .attr("style", "margin-left: auto; margin-right: auto");
  const xScale = d3
    .scaleLinear()
    .domain([Math.min(...x) / kx, Math.max(...x) / kx])
    .range([leftPadding, width - rightPadding]);
  const yScale = d3
    .scaleLinear()
    .domain([
      Math.min(...y.map((e) => Math.min(...e.filter((ee) => !isNaN(ee))))) / ky,
      Math.max(...y.map((e) => Math.max(...e.filter((ee) => !isNaN(ee))))) / ky,
    ])
    .range([height - bottomPadding, topPadding]);

  let xAxis = d3.axisBottom(xScale).tickSize(bottomPadding - height + topPadding);
  svg
    .append("g")
    .attr("class", "x-axis")
    .attr("transform", `translate(0, ${height - bottomPadding})`)
    .call(xAxis)
    .selectAll(".tick line")
    .attr("color", "lightgray")
    .select(".domain")
    .attr("stroke-width", "2");
  let yAxis = d3.axisLeft(yScale).tickSize(leftPadding - width + rightPadding);
  svg
    .append("g")
    .attr("class", "y-axis")
    .attr("transform", `translate(${leftPadding})`)
    .call(yAxis)
    .selectAll(".tick line")
    .attr("color", "lightgray")
    .select(".domain")
    .attr("stroke-width", "2");

  const yLine = d3
    .line()
    .x((d) => xScale(d.x))
    .y((d) => yScale(d.y))
    .curve(d3.curveLinear);
  for (let [i, yy] of y.entries()) {
    let yyd = yy.map((e, k) => ({ x: x[k] / kx, y: e / ky })).filter((e) => !isNaN(e.y));
    svg
      .append("path")
      .classed("data", true)
      .attr("d", yLine(yyd))
      .style("fill", "none")
      .style("stroke", colors[i])
      .style("stroke-width", "2");
  }

  svg
    .append("foreignObject")
    .attr("class", "x axis-label")
    .attr("width", width)
    .attr("height", bottomPadding)
    .attr("transform", `translate(${0}, ${height - bottomPadding / 2})`)
    .text("$L_{\\rm r}$ / μH");
  svg
    .append("foreignObject")
    .attr("class", "y axis-label")
    .attr("width", height)
    .attr("height", leftPadding)
    .attr(
      "transform",
      ` rotate(-90, ${height / 2}, ${leftPadding / 2})` +
        `translate(${-height / 2 + leftPadding / 2}, ${leftPadding / 2 - height / 2})`
    )
    .text(() => `${ylabel}`);
  MathJax.Hub.Queue(["Typeset", MathJax.Hub]);
}
