<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Run two File-coupled instances of the DummySolver, exchanging fields</description>
    <executable>foo</executable>
    <parameters>bar</parameters>
    <command>../DummySolver-g --verbose Dummy_3DCubeFile_0.xml cube.left.xml 1>>output.out 2>>output.err  & ../DummySolver-g --verbose Dummy_3DCubeFile_1.xml cube.right.xml 1>>output.out 2>>output.err &</command>
    <files>
        <file description="Session File 0">Dummy_3DCubeFile_0.xml</file>
        <file description="Session File 1">Dummy_3DCubeFile_1.xml</file>
        <file description="Mesh File">cube.left.xml</file>
        <file description="Mesh File">cube.right.xml</file>
        <file description="PETSc config file">.petscrc</file>
    </files>
    <metrics>
        <metric type="regex" id="0">
            <regex>^L 2 error\s*\(variable\s*\w0\w\)\s*:\s*([+-]?\d.+\d|-?\d|[+-]?nan|[+-]?inf).*</regex>
            <matches>
                <match>
                    <field variable="u0S" tolerance="1E-5">0.000122474</field>
                    <field variable="v0S" tolerance="1E-5">0.000122474</field>
                    <field variable="u0R" tolerance="1E-7">8.03482e-06</field>
                    <field variable="v0R" tolerance="1E-7">8.03482e-06</field>
                </match>
            </matches>
        </metric>
        <metric type="regex" id="1">
            <regex>^L 2 error\s*\(variable\s*\w1\w\)\s*:\s*([+-]?\d.+\d|-?\d|[+-]?nan|[+-]?inf).*</regex>
            <matches>
                <match>
                    <field variable="u1S" tolerance="1E-6">1.22474e-05</field>
                    <field variable="v1S" tolerance="1E-6">1.22474e-05</field>
                    <field variable="u1R" tolerance="1E-6">8.03713e-05</field>
                    <field variable="v1R" tolerance="1E-6">8.03713e-05</field>
                </match>
            </matches>
        </metric>
    </metrics>
</test>

