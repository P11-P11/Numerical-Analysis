using Images  
using FFTW
using StatsBase

"""
Como primer paso en el pipeline tenemos que pre procesar la imagen
Agregamos bordes con pixeles negros
"""

function rellenarNegro(imagen)
    n, m = size(imagen)
	N , M = n%16==0 ? n : n + (16 - n%16), m%16==0 ? m : m + (16 - m%16)

	matriz_negra = fill(RGB(0., 0., 0.), N, M)

	i ,j = (N - n)÷2 + 1, (M - m)÷2 + 1
	matriz_negra[i: i + n-1, j: j + m-1] = imagen

	return  matriz_negra
end

"
promediar vecinos son funciones auxiliares para pasar de 
RGB a Ycbcr y viceversa
Promedian/expanden de a grupos de cuatro
"

function promediarVecinos(mat)
	n,m = size(mat)
	mat2 = Array{Float64}(undef, n÷2, m÷2)
	
	for i in 1:n÷2
		for j in 1:m÷2
			casillero = 
				mat[2i-1,2j-1] + mat[2i,2j-1] + mat[2i-1,2j] + mat[2i,2j]
			
			mat2[i,j] = casillero/4
		end
	end
	
	return mat2
end


function expandirVecinos(mat)
	n,m = size(mat)
	mat2 = Array{Float64}(undef, n*2, m*2)
	for i in 1:2*n
		for j in 1:2*m
			
			mat2[i,j] = mat[i%2 + i÷2, j%2 + j÷2]
		end
	end
	return mat2
end

"
usamos el metodo YCBR de la libreria
"

function RGB_a_YCbCr(imagen)
	im = rellenarNegro(imagen)
	im1 = YCbCr.(im)
	
	channels = channelview(im1)
	Y  = channels[1,:,:]
	Cb = channels[2,:,:]
	Cr = channels[3,:,:]
	
	Cb = promediarVecinos(Cb)
	Cr = promediarVecinos(Cr)
	
	Y  = Y.-128
	Cb = Cb.-128
	Cr = Cr.-128
	
	return Y,Cb,Cr
end

function YCbCr_a_RGB(Y,cb,cr)
	Y   = Y.+128.
	cb  = cb.+128.
	cr  = cr.+128.

	return RGB.(colorview(YCbCr,Y,expandirVecinos(cb),expandirVecinos(cr)))
end

"
transforma coseno discreto de a bloques de cuatro 
"
function transformarBloques(mat)
	n,m = size(mat)
	
	for i in 1:8:n
		for j in 1:8:m
			dct!(@view(mat[i:i+7, j:j+7]))
		end
	end
	
	return mat
end
" idem pero inversa"
function transformarBloques⁻¹(mat)
	n,m = size(mat)
	for i in 1:8:n
		for j in 1:8:m
			idct!(@view(mat[i:i+7, j:j+7]))
		end
	end
	
	return mat
end
"
Quantizar toma bloques de a ocho y divide coordenada a coordenada
mij/Qij
"

function Quantizar(mat, Q)
	n, m = size(mat)
	res = Array{Float64}(undef, n, m)
	for i in (1:8:n) 
		for j in (1:8:m) 
			res[i:i+7,j:j+7] = round.(mat[i:i+7,j:j+7]./Q)
		end
	end
	return res
end

function Quantizar⁻¹(mat, Q)
	n, m = size(mat)
	res = Array{Float64}(undef, n, m)
	for i in (1:8:n) 
		for j in (1:8:m) 
			res[i:i+7,j:j+7] = (mat[i:i+7,j:j+7].*Q)
		end
	end
	return res
end

"zigzag es medio galerazo

lo entendemos con este ejemplo
si tenemos una matriz de 4*4

a11 a12 a13 a14
a21 a22 a23 a24
a31 a32 a33 a34
a41 a42 a43 a44

obersvar que para que elemento del zigzag los indices suman lo mismo
a11 // 1 + 1 = 2
a21 a12 // 2 + 1 = 1 + 2 = 3
a13 a22 a31 // 1 + 3 = 2 + 2 = 3 +1
a14 a23 a32 a41 // 1 + 4 = 2 + 3 = 3 + 2 = 4 + 1
...
entonces los elementos del zigzag quedan definidos por la suma de sus indices

"
function recorrido_zigzag(matriz)
    m, n = size(matriz)
    res = []
    
    for suma in 2:(m + n)
        if suma % 2 == 0
            for i in max(1, suma - m):min(suma - 1, n)
               push!(res, matriz[suma - i, i])

            end
		end
		if suma % 2 == 1
            for i in max(1, suma - n):min(suma - 1, m)
                push!(res,matriz[i, suma - i])

            end
        end
    end
    
    return res
end

function zigzag⁻¹(vector)
	n = 8
	m = 8
	bloque = fill(0, n, m)
    j = 1
    for suma in 2:(m + n)
        if suma % 2 == 1  
            for i in max(1, suma - n):min(suma - 1, m)
                bloque[i, suma - i] = vector[j]
                j += 1
            end
        else
            for i in max(1, suma - m):min(suma - 1, n)
                bloque[suma - i, i] = vector[j]
                j += 1
            end
        end
    end
    
    return bloque
end

function tiraRle(mat)
	n,m = size(mat)
	res = []
	for i in 1:8:n
		for j in 1:8:m
			vals, reps = rle(recorrido_zigzag(mat[i:i+7, j:j+7]))
			push!(res,reps...)
			push!(res,vals...)
		end
	end
	return res
end

function tiraRle⁻¹(rle,n,m,k)
	
	M = Array{Float64}(undef, n, m)
	a=1
	b=1
	
	N = size(rle)[1]
	sum = 0
	j = k
	i = k
	while((i <= N) && (a < n+1 && b<  m+1))
		sum += rle[i]
		if sum == 64
			reps = rle[j:i]
			vals = rle[i+1:i+1+i-j]
			
			vector = inverse_rle(vals,Int.(reps))
			bloque = zigzag⁻¹(vector)

			M[a:a+7,b:b+7] = bloque
			
			b=b+8
			if b>m
				b=1
				a=a+8
			end
			
			j = i+1 + i-j + 1 #i + largo intervalo +2 -> para el proximo rep
			i = j-1
			
			sum = 0
		end
		i += 1
	end
	return M, i
end

function de_imagen_a_rle(imagen,Q)#genera el vector rle a partir de la imagen y la matriz de cuantizacion
	Y,cb,cr = RGB_a_YCbCr(imagen)

	Y  = transformarBloques(Y)
	cb = transformarBloques(cb)
	cr = transformarBloques(cr)

	Y  = Quantizar(Y, Q)
	cb = Quantizar(cb, Q)
	cr = Quantizar(cr, Q) 
	
	return vcat(tiraRle(Y),tiraRle(cb),tiraRle(cr))
end


function de_rle_a_imagen(n,m,Q,rle)#reconstruye la imagen a partir del vector rle, la matriz de cuantizacion y el tamaño de la imagen
	
	Y , i_1 = tiraRle⁻¹(rle,n,m,1)
	cb, i_2 = tiraRle⁻¹(rle,n÷2,m÷2,i_1)
	cr, i_3 = tiraRle⁻¹(rle,n÷2, m÷2, i_2)
	
	Y  = Quantizar⁻¹(Y, Q)
	cb = Quantizar⁻¹(cb, Q)
	cr = Quantizar⁻¹(cr, Q)

	Y  = transformarBloques⁻¹(Y)
	cb = transformarBloques⁻¹(cb)
	cr = transformarBloques⁻¹(cr)

	return YCbCr_a_RGB(Y,cb,cr) 
end


function guardarDatos(name,n,m,Q,rle)

	io = open(name*".uba","w")

	write(io, UInt16[n])
	write(io, UInt16[m])
	
	#guardo la matriz Q
	for i in 1:8
		for j in 1:8
			write(io, UInt8[Q[i,j]])
		end
	end
	
	#guardo el vector rle
	for j in 1:length(rle)
		write(io, Int8[rle[j]])
	end
	
	close(io)
end

function leerDatos(name)
	io = open(name)

	n = read(io,UInt16)
	m = read(io,UInt16)
	
	Q   = Array{Float64}(undef, 8, 8)
	#leo matriz de cuantizacion
	for i in 1:8
		for j in 1:8
			Q[i,j] = read(io, UInt8)
		end
	end
	
	rle = []
	#leo todo el vector rle
	while !eof(io)
		push!(rle, read(io, Int8))
	end

	close(io)
	return n,m,Q,rle
end


function comprimir(name)#recibe el nombre de la imagen y guarda la compresion en un archivo de extension .uba
	imagen = load(name)

	#matriz de cuantizacion
	Q=[16 11 10 16 24 40 51 61; 
		12 12 14 19 26 58 60 55;
		14 13 16 24 40 57 69 56;
		14 17 22 29 51 87 80 62;
		18 22 37 56 68 109 103 77;
		24 35 55 64 81 104 113 92;
		49 64 78 87 103 121 120 101;
		72 92 95 98 112 100 103 99]
	
	n, m = size(imagen)
	
	
	guardarDatos(name,n,m,Q,de_imagen_a_rle(imagen,Q))
end

function descomprimir(archivo_uba)#recibe el archivo comprimido y devuelve la imagen
	
	n,m,Q,rle = leerDatos(archivo_uba)

	N , M = n%16==0 ? n : n + (16 - n%16), m%16==0 ? m : m + (16 - m%16)#CALCULO LOS N MOD 16
	
	imagen_con_bordes = de_rle_a_imagen(N,M,Q,rle)
	
	i ,j = (N - n)÷2 + 1, (M - m)÷2 + 1 #CACLULO LOS INDICES PARA LA MATRIZ SIN BORDES

	imagen_sin_borde= imagen_con_bordes[i: i + n-1, j: j + m-1]
	return imagen_sin_borde
end

# un test

#comprimir("paisaje.bmp")
#descomprimir("paisaje.bmp"*".uba")

function test()
	comprimir("tanque.jpg")
	img = descomprimir("tanque.jpg"*".uba")	
	save("compreso.bmp",img)
end
