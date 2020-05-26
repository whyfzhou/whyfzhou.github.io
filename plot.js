function drawSingleOP(data) {
  let oldPlot = document.querySelector("#plot");
  if (oldPlot) {
    oldPlot.remove();
  }

  let t = data.t;
  let ir = data.ir;
  let vr = data.vr;
  let im = data.im;
  let vhb = data.vhb;

  const padding = 30;
  const width = 400;
  const height = 300;

  const svg = d3
    .select("#chart")
    .append("svg")
    .attr("id", "plot")
    .attr("width", width)
    .attr("height", height);
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
  svg
    .append("g")
    .attr("class", "i-axis")
    .attr("transform", `translate(${padding})`)
    .call(iAxis);
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
}
