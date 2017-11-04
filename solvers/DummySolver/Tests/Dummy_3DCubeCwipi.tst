<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Run two CWIPI-coupled instances of the DummySolver, exchanging fields</description>
    <executable>foo</executable>
    <parameters>bar</parameters>
    <command>mpirun -n 1 ../DummySolver-g --verbose --cwipi 'Dummy0' Dummy_3DCubeCwipi_0.xml mesh.cube.xml : -n 1 ../DummySolver-g --verbose --cwipi 'Dummy1' Dummy_3DCubeCwipi_1.xml mesh.cube.xml 1>output.out 2>output.err</command>
    <files>
        <file description="Session File 0">Dummy_3DCubeCwipi_0.xml</file>
        <file description="Session File 1">Dummy_3DCubeCwipi_1.xml</file>
        <file description="Mesh File">mesh.cube.xml</file>
        <file description="PETSc config file">.petscrc</file>
    </files>
    <metrics>
        <metric type="regex" id="0">
            <regex>^L 2 error\s*\(variable\s*\w0\w\)\s*:\s*([+-]?\d.+\d|-?\d|[+-]?nan|[+-]?inf).*</regex>
            <matches>
                <match>
                    <field variable="u0S" tolerance="1E-5">0.00011547</field>
                    <field variable="v0S" tolerance="1E-5">0.00011547</field>
                    <field variable="u0R" tolerance="1E-6">1.03161e-05</field>
                    <field variable="v0R" tolerance="1E-6">1.03161e-05</field>
                </match>
            </matches>
        </metric>
        <metric type="regex" id="1">
            <regex>^L 2 error\s*\(variable\s*\w1\w\)\s*:\s*([+-]?\d.+\d|-?\d|[+-]?nan|[+-]?inf).*</regex>
            <matches>
                <match>
                    <field variable="u1S" tolerance="1E-6">1.1547e-05</field>
                    <field variable="v1S" tolerance="1E-6">1.1547e-05</field>
                    <field variable="u1R" tolerance="1E-5">0.000103161</field>
                    <field variable="v1R" tolerance="1E-5">0.000103161</field>
                </match>
            </matches>
        </metric>
    </metrics>
</test>

