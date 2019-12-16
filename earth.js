
//纹理
function handleLoadedTexture(texture) {
  gl.pixelStorei(gl.UNPACK_FLIP_Y_WEBGL, true);
  gl.bindTexture(gl.TEXTURE_2D, texture);
  gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, gl.RGBA, gl.UNSIGNED_BYTE, texture.image);
  gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.LINEAR);
  gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.LINEAR_MIPMAP_NEAREST);
  gl.generateMipmap(gl.TEXTURE_2D);

  gl.bindTexture(gl.TEXTURE_2D, null);
}
var earthTexture;
function initTexture(src) {
  earthTexture = gl.createTexture();
  earthTexture.image = new Image();
  earthTexture.image.onload = function () {
    handleLoadedTexture(earthTexture)
  }

  earthTexture.image.src = src;
}

// function mvPushMatrix() {
//   var copy = ""//new okMat4(); //TODO
//   mvMatrix.clone(copy);
//   mvMatrixStack.push(copy);
// }
// function mvPopMatrix() {
//   if (mvMatrixStack.length == 0) {
//     throw "Invalid popMatrix!";
//   }
//   mvMatrix = mvMatrixStack.pop();
// }

// function setMatrixUniforms() {
//   gl.uniformMatrix4fv(gl.program.pMatrixUniform, false, pMatrix);
//   gl.uniformMatrix4fv(gl.program.mvMatrixUniform, false, mvMatrix);

//   //var normalMatrix = mvMatrix.inverse().transpose();
//   var normalMatrix = Mat4.identity(Mat4.create());
//   Mat4.inverse(normalMatrix, mvMatrix);
//   Mat4.transpose(normalMatrix, normalMatrix);

//   gl.uniformMatrix4fv(gl.program.nMatrixUniform, false, normalMatrix);
// }





//var earthVertexPositionBuffer;
var earthVertexNormalBuffer;
var earthVertexTextureCoordBuffer;
var earthVertexIndexBuffer;
//球体的三角剖分 http://lulier.github.io/2018/06/12/webGL%E7%90%83%E4%BD%93%E6%A8%A1%E5%9E%8B/index.html
function initBuffers() {
  var latitudeBands = 30;
  var longitudeBands = 30;
  var radius = 2;

  var vertexPositionData = [];
  var normalData = [];
  var textureCoordData = [];
  for (var latNumber = 0; latNumber <= latitudeBands; latNumber++) {
    var theta = latNumber * Math.PI / latitudeBands;
    var sinTheta = Math.sin(theta);
    var cosTheta = Math.cos(theta);

    for (var longNumber = 0; longNumber <= longitudeBands; longNumber++) {
      var phi = longNumber * 2 * Math.PI / longitudeBands;
      var sinPhi = Math.sin(phi);
      var cosPhi = Math.cos(phi);

      var x = cosPhi * sinTheta;
      var y = cosTheta;
      var z = sinPhi * sinTheta;
      var u = 1 - (longNumber / longitudeBands);
      var v = 1 - (latNumber / latitudeBands);

      normalData.push(x);
      normalData.push(y);
      normalData.push(z);
      textureCoordData.push(u);
      textureCoordData.push(v);
      vertexPositionData.push(radius * x);
      vertexPositionData.push(radius * y);
      vertexPositionData.push(radius * z);
    }
  }

  var indexData = [];
  for (var latNumber = 0; latNumber < latitudeBands; latNumber++) {
    for (var longNumber = 0; longNumber < longitudeBands; longNumber++) {
      var first = (latNumber * (longitudeBands + 1)) + longNumber;
      var second = first + longitudeBands + 1;
      indexData.push(first);
      indexData.push(second);
      indexData.push(first + 1);

      indexData.push(second);
      indexData.push(second + 1);
      indexData.push(first + 1);
    }
  }

  earthVertexNormalBuffer = gl.createBuffer();
  gl.bindBuffer(gl.ARRAY_BUFFER, earthVertexNormalBuffer);
  gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(normalData), gl.STATIC_DRAW);
  earthVertexNormalBuffer.itemSize = 3;
  earthVertexNormalBuffer.numItems = normalData.length / 3;

  earthVertexTextureCoordBuffer = gl.createBuffer();
  gl.bindBuffer(gl.ARRAY_BUFFER, earthVertexTextureCoordBuffer);
  gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(textureCoordData), gl.STATIC_DRAW);
  earthVertexTextureCoordBuffer.itemSize = 2;
  earthVertexTextureCoordBuffer.numItems = textureCoordData.length / 2;

  earthVertexPositionBuffer = gl.createBuffer();
  gl.bindBuffer(gl.ARRAY_BUFFER, earthVertexPositionBuffer);
  gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(vertexPositionData), gl.STATIC_DRAW);
  earthVertexPositionBuffer.itemSize = 3;
  earthVertexPositionBuffer.numItems = vertexPositionData.length / 3;
  //initArrayBuffer(gl, 'aVertexPosition', vertexPositionData, 3, gl.FLOAT);

  // earthVertexIndexBuffer = gl.createBuffer();
  // gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, earthVertexIndexBuffer);
  // gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, new Uint16Array(indexData), gl.STATIC_DRAW);
  // earthVertexIndexBuffer.itemSize = 1;
  // earthVertexIndexBuffer.numItems = indexData.length;
  IDXS_LEN = indexData.length;
  initIBOBuffer(gl, indexData, gl.UNSIGNED_SHORT);

  return {
    i: indexData
  }
}


function drawScene() {
  gl.viewport(0, 0, gl.viewportWidth, gl.viewportHeight);
  gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);

  //pMatrix = okMat4Proj(45, gl.viewportWidth / gl.viewportHeight, 0.1, 100.0);//TODO
  Mat4.perspective(pMatrix, 45, gl.viewportWidth / gl.viewportHeight, 0.1, 100);

  var lighting = 0;
  gl.uniform1i(gl.program.useLightingUniform, lighting);
  // if (lighting) {
  //     gl.uniform3f(
  //         gl.program.ambientColorUniform,
  //         parseFloat(document.getElementById("ambientR").value),
  //         parseFloat(document.getElementById("ambientG").value),
  //         parseFloat(document.getElementById("ambientB").value)
  //     );

  //     var lightingDirection = new okVec3(
  //         parseFloat(document.getElementById("lightDirectionX").value),
  //         parseFloat(document.getElementById("lightDirectionY").value),
  //         parseFloat(document.getElementById("lightDirectionZ").value)
  //     );

  //     var adjustedLD = lightingDirection.normalize(false);
  //     adjustedLD = okVec3MulVal(adjustedLD, -1.0);
  //     gl.uniform3fv(gl.program.lightingDirectionUniform, adjustedLD.toArray());

  //     gl.uniform3f(
  //         gl.program.directionalColorUniform,
  //         parseFloat(document.getElementById("directionalR").value),
  //         parseFloat(document.getElementById("directionalG").value),
  //         parseFloat(document.getElementById("directionalB").value)
  //     );
  // }

  //mvMatrix = okMat4Trans(0.0, 0.0, -6.0); //TODO
  //mvMatrix = okMat4Mul(mvMatrix, earthRotationMatrix);//TODO
  //Mat4.translate(mvMatrix, mvMatrix, [0.0, 0.0, -6.0]);
  var look = 5;
  Mat4.lookAt(mvMatrix, [0, look, look], [0, 0, 0], [0, 1, 0]);
  Mat4.multiply(mvMatrix, mvMatrix, earthRotationMatrix);

  gl.activeTexture(gl.TEXTURE0);
  gl.bindTexture(gl.TEXTURE_2D, earthTexture);
  gl.uniform1i(gl.program.samplerUniform, 0);

  gl.bindBuffer(gl.ARRAY_BUFFER, earthVertexPositionBuffer);
  gl.vertexAttribPointer(gl.program.vertexPositionAttribute, earthVertexPositionBuffer.itemSize, gl.FLOAT, false, 0, 0);

  gl.bindBuffer(gl.ARRAY_BUFFER, earthVertexTextureCoordBuffer);
  gl.vertexAttribPointer(gl.program.textureCoordAttribute, earthVertexTextureCoordBuffer.itemSize, gl.FLOAT, false, 0, 0);

  gl.bindBuffer(gl.ARRAY_BUFFER, earthVertexNormalBuffer);
  gl.vertexAttribPointer(gl.program.vertexNormalAttribute, earthVertexNormalBuffer.itemSize, gl.FLOAT, false, 0, 0);

  //gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, earthVertexIndexBuffer);


  //setMatrixUniforms();
  gl.uniformMatrix4fv(gl.program.pMatrixUniform, false, pMatrix);
  gl.uniformMatrix4fv(gl.program.mvMatrixUniform, false, mvMatrix);
  var normalMatrix = Mat4.identity(Mat4.create());
  Mat4.inverse(normalMatrix, mvMatrix);
  Mat4.transpose(normalMatrix, normalMatrix);
  gl.uniformMatrix4fv(gl.program.nMatrixUniform, false, normalMatrix);

  gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
  gl.drawElements(gl.TRIANGLES, IDXS_LEN, gl.UNSIGNED_SHORT, 0);
}


function render() {
  requestAnimationFrame(render);
  drawScene();
}


var gl;
var pMatrix;
var mvMatrix;
var mvMatrixStack = [];
var IDXS_LEN = 0;
Mat4 = new mat4();
function loadWebGL() {
  gl = getWebGLContext("canvasId");
  var vertexShader = loadShader(gl, gl.VERTEX_SHADER, VSHADER_SOURCE);
  var fragmentShader = loadShader(gl, gl.FRAGMENT_SHADER, FSHADER_SOURCE);
  gl.program = createProgram(vertexShader, fragmentShader);

  //vbo
  gl.program.vertexPositionAttribute = gl.getAttribLocation(gl.program, "aVertexPosition");
  gl.enableVertexAttribArray(gl.program.vertexPositionAttribute);

  gl.program.textureCoordAttribute = gl.getAttribLocation(gl.program, "aTextureCoord");
  gl.enableVertexAttribArray(gl.program.textureCoordAttribute);

  gl.program.vertexNormalAttribute = gl.getAttribLocation(gl.program, "aVertexNormal");
  gl.enableVertexAttribArray(gl.program.vertexNormalAttribute);

  //获取变量
  gl.program.pMatrixUniform = gl.getUniformLocation(gl.program, "uPMatrix");
  gl.program.mvMatrixUniform = gl.getUniformLocation(gl.program, "uMVMatrix");
  gl.program.nMatrixUniform = gl.getUniformLocation(gl.program, "uNMatrix");
  gl.program.samplerUniform = gl.getUniformLocation(gl.program, "uSampler");
  gl.program.useLightingUniform = gl.getUniformLocation(gl.program, "uUseLighting");
  gl.program.ambientColorUniform = gl.getUniformLocation(gl.program, "uAmbientColor");
  gl.program.lightingDirectionUniform = gl.getUniformLocation(gl.program, "uLightingDirection");
  gl.program.directionalColorUniform = gl.getUniformLocation(gl.program, "uDirectionalColor");

  //初始mat
  pMatrix = Mat4.identity(Mat4.create());
  mvMatrix = Mat4.identity(Mat4.create());
  earthRotationMatrix = Mat4.identity(Mat4.create());


  initBuffers();
  initTexture("./world.png");


  gl.clearColor(0.0, 0.0, 0.0, 1.0);
  gl.enable(gl.DEPTH_TEST);
  //gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
  render();

  //初始化交互
  initInteracting();
}



function initInteracting() {
  var canvas = document.getElementById("canvasId");
  canvas.onmousedown = handleMouseDown;
  document.onmouseup = handleMouseUp;
  document.onmousemove = handleMouseMove;
}
/////////////////////////////////////////////////
//      交互    
/////////////////////////////////////////////////
var mouseDown = false;
var lastMouseX = null;
var lastMouseY = null;

var orgMouseX = null;
var orgMouseY = null;
//var earthRotationMatrix = new Matrix4();
function handleMouseDown(event) {
  mouseDown = true;
  lastMouseX = event.clientX;
  lastMouseY = event.clientY;
  if (!orgMouseX) {
    orgMouseX = lastMouseX;
    orgMouseY = lastMouseY;
  }
}
function handleMouseUp(event) {
  mouseDown = false;
}

function handleMouseMove(event) {
  if (!mouseDown) {
    return;
  }
  var newX = event.clientX;
  var newY = event.clientY;
  var speedFactor = 60;

  var deltaX = newX - lastMouseX
  var newRotationMatrix = Mat4.identity(Mat4.create());//new okMat4();//TODO
  Mat4.rotate(newRotationMatrix, newRotationMatrix, (deltaX / speedFactor), [0, 1, 0]);//Y axis

  var deltaY = newY - lastMouseY;
  Mat4.rotate(newRotationMatrix, newRotationMatrix, (deltaY / speedFactor), [1, 0, 0]);//X axis

  Mat4.multiply(earthRotationMatrix, newRotationMatrix, earthRotationMatrix);

  lastMouseX = newX
  lastMouseY = newY;
}
/////////////////////////////////////////////////