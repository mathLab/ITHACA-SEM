<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Gmsh high-order tri square, order 8</description>
    <executable>MeshConvert</executable>
    <parameters>-m jac:list SquareTri.msh SquareTri.xml:xml:test</parameters>
    <files>
        <file description="Input File">SquareTri.msh</file>
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
