<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Run a single instances of the DummySolver, coupled with itself</description>
    <executable>DummySolver</executable>
    <parameters>Dummy_3DCubeFileSingle.xml cube.xml</parameters>
    <files>
        <file description="Session File">Dummy_3DCubeFileSingle.xml</file>
        <file description="Mesh File">cube.xml</file>
        <file description="PETSc config file">.petscrc</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u0S" tolerance="1e-6">4.89898e-05</value>
            <value variable="v0S" tolerance="1e-6">4.89898e-05</value>
            <value variable="u0R" tolerance="1e-6">3.26599e-05</value>
            <value variable="v0R" tolerance="1e-6">3.26599e-05</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u0S" tolerance="1e-6">9e-05</value>
            <value variable="v0S" tolerance="1e-6">9e-05</value>
            <value variable="u0R" tolerance="1e-6">6e-05</value>
            <value variable="v0R" tolerance="1e-6">6e-05</value>
        </metric>
    </metrics>
</test>
