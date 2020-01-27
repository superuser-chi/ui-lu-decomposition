module.exports = {
  pluginOptions: {
    quasar: {
      importStrategy: "kebab",
      rtlSupport: false
    }
  },
  transpileDependencies: ["quasar"],
  devServer: {
    proxy: {
      "/api": {
        target: "http://pdftexapi.pythonanywhere.com",
        changeOrigin: true,
        ws: true,
        pathRewrite: {
          "^/api": ""
        }
      },
      "/faas": {
        target: "https://faas.serv.mathsgaze.com/function/",
        changeOrigin: true,
        ws: true,
        pathRewrite: {
          "^/faas": ""
        }
      }
    }
  }
};
