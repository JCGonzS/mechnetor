
function makeProtBarplot(id, proteins, interactors, max) {
  var my_chart = echarts.init(document.getElementById(id));
  var int_per = []
  interactors.forEach(function(i) {
    int_per.push(i/max*100);
  });
  // MAKE IT VARIABLE AND NETWORK-DEPENDANT
  // var nodes = cy.$("node[role='protein_main']");
  // var max_ints = nodes.length - 1;
  // var proteins = [];
  // var ints = [];
  // nodes.forEach(function(node) {
  //   proteins.push(node.data("label"));
  // });
  var option = {
      title: {
          text: 'Proteins ranked by number of interactors',
          subtext: "Maximum number of interactors per protein: "+max
      },
      grid: {
        containLabel: true,
        bottom: "5%"
      },
      tooltip: {
        trigger: 'axis',
        axisPointer: {
          type: 'shadow'
        }
      },
      xAxis: {
        minInterval: 1
      },
      yAxis: {
          type: 'category',
          data: proteins,
          axisLabel: {fontWeight: 'bold'}
      },
      series: [{
          name: 'Interactors',
          type: 'bar',
          color: ['#5d6d7e'],
          data: interactors
      }]
    };
    my_chart.setOption(option);
  };

function makeIntTypeBarplot(id, labels, series, max) {
  var my_chart = echarts.init(document.getElementById(id));

  var option = {
      title: {
          text: 'Interaction types ranked',
          subtext: "by number of linked protein pairs. Max: "+max
      },
      grid: {
        containLabel: true,
        bottom: "5%"
      },
      tooltip: {
        trigger: 'axis',
        axisPointer: {
          type: 'shadow'
        }
      },
      xAxis: {
        minInterval: 1
      },
      yAxis: {
          type: 'category',
          data: labels,
          axisLabel: {fontWeight: 'bold'}
      },
      series: [{
          name: 'Linked protein pairs',
          type: 'bar',
          data: series
      }]
  };
  my_chart.setOption(option);
};

function makeDomBarplot(id, doms, dom_ints_tot, max) {
  var my_chart = echarts.init(document.getElementById(id));
  var option = {
      title: {
          text: 'Domains ranked by number of interactors',
          subtext: "Maximum number of interactors per protein: "+max,
          textAlign: "center"
      },
      grid: {
        left: '15%',
        bottom: "5%",
        containLabel: true
      },
      tooltip: {
        trigger: 'axis',
        axisPointer: {
          type: 'shadow'
        }
      },
      xAxis: {
        minInterval: 1
      },
      yAxis: {
          type: 'category',
          data: doms,
          axisLabel: {fontWeight: 'bold'}
      },
      series: [{
          name: 'Total interactors',
          type: 'bar',
          color: ['#717d7e'],
          data: dom_ints_tot
      }]
    };
    my_chart.setOption(option);
};

function makeDomBarplot2(id, doms, dom_ints, max) {
  var my_chart = echarts.init(document.getElementById(id));
  var option = {
      title: {
          text: 'Domains ranked by number of interactors',
          subtext: "Maximum number of interactors per protein: "+max,
          textAlign: "center"
      },
      legend: {
        data: ["DDI", "iDDI", "DMI"]
      },
      grid: {
        left: '15%',
        bottom: "5%",
        containLabel: true
      },
      tooltip: {
        trigger: 'axis',
        axisPointer: {
          type: 'shadow'
        }
      },
      xAxis: {
        minInterval: 1
      },
      yAxis: {
          type: 'category',
          data: doms,
          axisLabel: {fontWeight: 'bold'}
      },
      series: [{
          name: 'DDI interactors',
          type: 'bar',
          stack: 'ints',
          // label: {
          //   show: true,
          //   position: "insideRight"
          // },
          color: ['#16A085'],
          data: dom_ints.DDI
      }, {
        name: 'iDDI interactors',
        type: 'bar',
        stack: 'ints',
        // label: {
        //   show: true,
        //   position: "insideRight"
        // },
        color: ['#D4AC0D'],
        data: dom_ints.iDDI
      }, {
        name: 'DMI interactors',
        type: 'bar',
        stack: 'ints',
        // label: {
        //   show: true,
        //   position: "insideRight"
        // },
        color: ['#AF7AC5'],
        data: dom_ints.DMI
      }]
    };
    my_chart.setOption(option);
};
