import struct
import math

def parse_flags(high_byte, low_byte):
    """
    Parse the high and low bytes for color and flags.
    """
    color1 = (high_byte & 0xF0) >> 4
    color2 = high_byte & 0x0F
    shape = 'Triangle' if (low_byte & 0x10) else 'Square'
    backface_culling = (low_byte & 0x40) != 0
    return color1, color2, shape, backface_culling

def read_vertex(file):
    """
    Reads a vertex consisting of three signed words (X, Y, Z) without endian-ness conversion.
    """
    data = file.read(6)
    if len(data) < 6:
        raise EOFError("Incomplete vertex data.")
    x, y, z = struct.unpack('>hhh', data)
    return (x, y, z)

def compute_normal(v1, v2, v3):
    """
    Compute the unit normal vector for a face given three vertices.
    """
    ux, uy, uz = v2[0] - v1[0], v2[1] - v1[1], v2[2] - v1[2]
    vx, vy, vz = v3[0] - v1[0], v3[1] - v1[1], v3[2] - v1[2]

    nx = uy * vz - uz * vy
    ny = uz * vx - ux * vz
    nz = ux * vy - uy * vx

    length = math.sqrt(nx * nx + ny * ny + nz * nz)
    if length == 0:
        return (0.0, 0.0, 0.0)
    return (nx / length, ny / length, nz / length)

def compute_2d_bounding_box(vertices):
    """
    Compute a 2D bounding box in the x/z plane for a given set of vertices.
    """
    xs = [v[0] for v in vertices]
    zs = [v[2] for v in vertices]
    return (min(xs), max(xs), min(zs), max(zs))

def convert_models_to_c_3d_arrays(input_file, c_output_file):
    """
    Reads as many well-formed models as possible from the input binary file and outputs them
    into 3D C arrays, including 3D arrays for face normals and a 2D bounding box for each model.
    """
    models_vertices = []
    models_faces = []
    models_colors = []
    models_normals = []
    models_bounding_boxes = []
    vertex_counts = []
    face_counts = []
    color_counts = []

    with open(input_file, 'rb') as f:
        while True:
            try:
                header = f.read(2)
                if len(header) < 2:
                    break

                poly_count_minus_one = struct.unpack('>H', header)[0]
                poly_count = poly_count_minus_one + 1

                vertices, faces, colors, normals = [], [], [], []

                for _ in range(poly_count):
                    poly_data = f.read(2)
                    if len(poly_data) < 2:
                        raise EOFError("Incomplete polygon data.")

                    high_byte, low_byte = struct.unpack('>BB', poly_data)
                    color1, color2, shape, backface_culling = parse_flags(high_byte, low_byte)
                    vertex_count = 3 if shape == 'Triangle' else 4
                    current_poly_vertices = []

                    for _ in range(vertex_count):
                        vertex = read_vertex(f)
                        if vertex not in vertices:
                            vertices.append(vertex)
                        current_poly_vertices.append(vertices.index(vertex))

                    if not backface_culling:
                        current_poly_vertices.reverse()

                    if vertex_count == 3:
                        faces.append(current_poly_vertices)
                        colors.append((color1, color2))
                        v1, v2, v3 = [vertices[idx] for idx in current_poly_vertices]
                        normals.append(compute_normal(v1, v2, v3))
                    else:
                        faces.append([current_poly_vertices[0], current_poly_vertices[1], current_poly_vertices[2]])
                        faces.append([current_poly_vertices[0], current_poly_vertices[2], current_poly_vertices[3]])
                        colors.extend([(color1, color2), (color1, color2)])
                        v1, v2, v3, v4 = [vertices[idx] for idx in current_poly_vertices]
                        normals.append(compute_normal(v1, v2, v3))
                        normals.append(compute_normal(v1, v3, v4))

                models_vertices.append(vertices)
                models_faces.append(faces)
                models_colors.append(colors)
                models_normals.append(normals)
                models_bounding_boxes.append(compute_2d_bounding_box(vertices))

                vertex_counts.append(len(vertices))
                face_counts.append(len(faces))
                color_counts.append(len(colors))

            except EOFError:
                print("Reached end of file or incomplete model data.")
                break

    with open(c_output_file, 'w') as out_c:
        out_c.write("// Generated 3D C arrays for models with face normals\n\n")

        out_c.write(f"int models_vertices[{len(models_vertices)}][{max(len(v) for v in models_vertices)}][3] = {{\n")
        for vertices in models_vertices:
            out_c.write("    {\n")
            for v in vertices:
                out_c.write(f"        {{{v[0]}, {v[1]}, {v[2]}}},\n")
            out_c.write("    },\n")
        out_c.write("};\n\n")

        out_c.write(f"int models_faces[{len(models_faces)}][{max(len(f) for f in models_faces)}][3] = {{\n")
        for faces in models_faces:
            out_c.write("    {\n")
            for face in faces:
                out_c.write(f"        {{{face[0]}, {face[1]}, {face[2]}}},\n")
            out_c.write("    },\n")
        out_c.write("};\n\n")

        out_c.write(f"int models_colors[{len(models_colors)}][{max(len(c) for c in models_colors)}][2] = {{\n")
        for colors in models_colors:
            out_c.write("    {\n")
            for color in colors:
                out_c.write(f"        {{{color[0]}, {color[1]}}},\n")
            out_c.write("    },\n")
        out_c.write("};\n\n")

        out_c.write(f"float models_normals[{len(models_normals)}][{max(len(n) for n in models_normals)}][3] = {{\n")
        for normals in models_normals:
            out_c.write("    {\n")
            for normal in normals:
                out_c.write(f"        {{{normal[0]:.6f}, {normal[1]:.6f}, {normal[2]:.6f}}},\n")
            out_c.write("    },\n")
        out_c.write("};\n\n")

        out_c.write(f"int vertex_counts[{len(vertex_counts)}] = {{\n")
        for count in vertex_counts:
            out_c.write(f"    {count},\n")
        out_c.write("};\n\n")

        out_c.write(f"int face_counts[{len(face_counts)}] = {{\n")
        for count in face_counts:
            out_c.write(f"    {count},\n")
        out_c.write("};\n\n")

        out_c.write(f"int models_bounding_boxes[{len(models_bounding_boxes)}][4] = {{\n")
        for bbox in models_bounding_boxes:
            out_c.write(f"    {{{bbox[0]}, {bbox[1]}, {bbox[2]}, {bbox[3]}}},\n")
        out_c.write("};\n\n")

    print(f"All models converted with face normals and bounding boxes to 3D C arrays in file: {c_output_file}.")

convert_models_to_c_3d_arrays('model3d.bin', 'models_3d_data.h')