<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Simple geometry with sphere</description>
    <executable>NekMesh</executable>
    <parameters>-m jac:list 3d_sphere.mcf 3d_sphere-out.xml:xml:test</parameters>
    <files>
        <file description="Input File">3d_sphere.mcf</file>
        <file description="Input File 2">3d_sphere.stp</file>
    </files>
    <metrics>
        <metric type="regex" id="1">
            <regex>^Total negative Jacobians: (\d+)</regex>
            <matches>
                <match>
                    <field id="0">0</field>
                </match>
            </matches>
        </metric>
    </metrics>
</test>
