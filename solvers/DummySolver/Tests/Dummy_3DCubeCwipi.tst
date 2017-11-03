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
        <metric type="L2" id="1">
            <value variable="aS" tolerance="1E-5">0.00011547</value>
            <value variable="bS" tolerance="1E-5">0.00011547</value>
            <value variable="aR" tolerance="1E-6">1.03161e-05</value>
            <value variable="bR" tolerance="1E-6">1.03161e-05</value>
        </metric>
        <metric type="L2" id="2">
            <value variable="cS" tolerance="1E-6">1.1547e-05</value>
            <value variable="dS" tolerance="1E-6">1.1547e-05</value>
            <value variable="cR" tolerance="1E-5">0.000103161</value>
            <value variable="dR" tolerance="1E-5">0.000103161</value>
        </metric>
    </metrics>
</test>
