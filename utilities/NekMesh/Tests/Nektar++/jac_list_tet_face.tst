<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Test with an invalid face curvature</description>
    <executable>NekMesh</executable>
    <parameters>-m jac:list jac_list_tet_face.xml jac_list_tet_face-out.xml:xml:test</parameters>
    <files>
        <file description="Input File">jac_list_tet_face.xml</file>
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
