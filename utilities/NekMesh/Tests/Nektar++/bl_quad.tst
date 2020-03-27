<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Test BL refinement in 2D from XML file</description>
    <executable>NekMesh</executable>
    <parameters>-m bl:surf=74 -m jac:list bl_quad.xml bl_quad-out.xml:xml:test</parameters>
    <files>
        <file description="Input File">bl_quad.xml</file>
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
